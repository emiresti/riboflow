/* Pyodide WebAssembly Bridge for GitHub Pages static deployment.
 * Intercepts /api/* fetch calls and routes them into the exact same Python
 * backend package running inside the visitor's browser via Pyodide (CPython 3.12).
 */
(function () {
  const PYODIDE_INDEX = "https://cdn.jsdelivr.net/pyodide/v0.27.5/full/";
  const APP_TITLE = "RiboFlow";
  const origFetch = window.fetch.bind(window);

  // Create sleek non-intrusive boot overlay while Pyodide loads
  const overlay = document.createElement("div");
  overlay.id = "pyodide-boot-overlay";
  overlay.innerHTML = `
    <div style="
      position: fixed; inset: 0; z-index: 99999;
      background: rgba(8, 13, 22, 0.92); backdrop-filter: blur(8px);
      display: flex; align-items: center; justify-content: center;
      font-family: system-ui, -apple-system, sans-serif; color: #e6eefb;
    ">
      <div style="
        max-width: 440px; width: 90%; padding: 26px 28px; border-radius: 16px;
        background: #0f1726; border: 1px solid #2dd4bf55;
        box-shadow: 0 24px 60px rgba(0, 0, 0, 0.55); text-align: center;
      ">
        <div style="font-size: 2rem; margin-bottom: 10px;">🧬</div>
        <h2 style="margin: 0 0 8px; font-size: 1.2rem; color: #2dd4bf;">
          ${APP_TITLE} — WebAssembly Python Lab
        </h2>
        <p id="pyodide-boot-msg" style="margin: 0 0 16px; font-size: 0.9rem; color: #94a3b8; line-height: 1.45;">
          Loading Python bioinformatics engine in your browser...<br/>
          <span style="font-size: 0.82rem; opacity: 0.85;">Tarayıcınızda Python biyoinformatik motoru başlatılıyor (~3-5 sn)...</span>
        </p>
        <div style="height: 6px; border-radius: 999px; background: #1e293b; overflow: hidden;">
          <div id="pyodide-boot-bar" style="
            width: 28%; height: 100%; border-radius: 999px;
            background: linear-gradient(90deg, #2dd4bf, #38bdf8);
            transition: width 0.35s ease;
          "></div>
        </div>
      </div>
    </div>
  `;

  function mountOverlay() {
    if (document.body && !document.getElementById("pyodide-boot-overlay")) {
      document.body.appendChild(overlay);
    }
  }
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", mountOverlay);
  } else {
    mountOverlay();
  }

  function setProgress(pct, msg) {
    const bar = document.getElementById("pyodide-boot-bar");
    const txt = document.getElementById("pyodide-boot-msg");
    if (bar) bar.style.width = pct + "%";
    if (txt && msg) txt.innerHTML = msg;
  }

  function loadScript(src) {
    return new Promise((resolve, reject) => {
      const s = document.createElement("script");
      s.src = src;
      s.onload = resolve;
      s.onerror = () => reject(new Error("Failed to load " + src));
      document.head.appendChild(s);
    });
  }

  async function bootPyodide() {
    mountOverlay();
    setProgress(35);
    const [bundleResp] = await Promise.all([
      origFetch("./py_bundle.json"),
      loadScript(PYODIDE_INDEX + "pyodide.js"),
    ]);
    if (!bundleResp.ok) {
      throw new Error("Could not load py_bundle.json (" + bundleResp.status + ")");
    }
    const bundle = await bundleResp.json();
    setProgress(65);

    const pyodide = await window.loadPyodide({ indexURL: PYODIDE_INDEX });
    setProgress(85);

    // Mount all Python source files into virtual filesystem
    for (const [relPath, content] of Object.entries(bundle.files)) {
      const full = "/app/src/" + relPath;
      const dir = full.substring(0, full.lastIndexOf("/"));
      pyodide.FS.mkdirTree(dir);
      pyodide.FS.writeFile(full, content, { encoding: "utf8" });
    }

    pyodide.runPython("import sys\nif '/app/src' not in sys.path:\n    sys.path.insert(0, '/app/src')\n");
    pyodide.runPython(bundle.adapter);
    setProgress(100);

    const el = document.getElementById("pyodide-boot-overlay");
    if (el) el.remove();
    return pyodide;
  }

  const pyodideReady = bootPyodide().catch((err) => {
    console.error("Pyodide boot error:", err);
    setProgress(100, `<span style="color:#f43f5e;">Engine failed to load: ${err.message}</span>`);
    throw err;
  });

  window.fetch = async function (input, init) {
    const rawUrl = typeof input === "string" ? input : (input && input.url) || "";
    let parsed;
    try {
      parsed = new URL(rawUrl, window.location.href);
    } catch {
      return origFetch(input, init);
    }

    const apiIdx = parsed.pathname.indexOf("/api/");
    if (apiIdx === -1) {
      return origFetch(input, init);
    }

    const apiPath = parsed.pathname.slice(apiIdx) + (parsed.search || "");
    const method = ((init && init.method) || (typeof input !== "string" && input.method) || "GET").toUpperCase();
    let bodyStr = "";
    if (init && init.body != null) {
      bodyStr = typeof init.body === "string" ? init.body : String(init.body);
    }

    const pyodide = await pyodideReady;
    const dispatchFn = pyodide.globals.get("_pyodide_dispatch");
    const resultJson = dispatchFn(method, apiPath, bodyStr);
    const resObj = JSON.parse(resultJson);

    return new Response(resObj.body, {
      status: resObj.status,
      headers: { "Content-Type": resObj.content_type || "application/json; charset=utf-8" },
    });
  };
})();
