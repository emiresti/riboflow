/* RiboFlow Interactive Web Laboratory */

const I18N = {
  tr: {
    subtitle: "Eş-Zamanlı Ribozom Trafiği (ℓ-TASEP) ve Kodon Harmonizasyonu Laboratuvarı",
    sampleLabel: "Hazır Biyolojik Gen / CDS",
    hostLabel: "Çeviri (Ekspresyon) Konak Organizması",
    alphaLabel: "Başlatma Hızı (α)",
    runSim: "Simülasyonu Çalıştır",
    customFasta: "Özel FASTA Yapıştır",
    pasteFasta: "Kodlayan DNA / mRNA dizisini yapıştırın (FASTA veya ham ATG...TAA):",
    loadSeq: "Diziyi Analiz Et ve Simüle Et",
    kpiFlux: "Protein Üretim Akışı (J)",
    kpiFluxSub: "protein / dakika",
    kpiCollision: "Disome Çarpışma Süresi",
    kpiNgd: "NGD / RQC Risk Skoru",
    kpiCai: "CAI ve GC / GC3",
    kpiPhase: "TASEP Trafik Fazı",
    tabTraffic: "1. Canlı Ribozom Trafiği & Ribo-seq",
    tabMutation: "2. Sessiz Mutasyon Laboratuvarı",
    tabHarmonize: "3. Kodon Harmonizatörü",
    stageTitle: "Gerçek Zamanlı ℓ-TASEP Ribozom Çeviri Görselleştiricisi",
    stageHint: "Her ribozom mRNA üzerinde 10 kodon (~30 nt) kaplar. Kırmızı halkalar No-Go Decay (NGD) tetikleyen disome çarpışmalarını gösterir.",
    speedLabel: "Hız",
    legFast: "Hızlı kodon",
    legMid: "Orta hızda kodon",
    legSlow: "Nadir / duraklama kodonu",
    legRibo: "Aktif ribozom",
    legColl: "Çarpışan disome (NGD riski)",
    chartRiboseqTitle: "Simüle Ribo-seq & Disome-seq Profili",
    chartRiboseqSub: "mRNA boyunca zaman-ortalamalı ribozom yoğunluğu (turkuaz) ve çarpışma tepeleri (kırmızı).",
    chartMinmaxTitle: "%MinMax Eş-Zamanlı Katlanma Ritmi",
    chartMinmaxSub: "Pozitif tepeler hızlı bölgeleri, negatif amber çukurları protein katlanma duraklamalarını gösterir.",
    mutTitle: 'Sinonim ("Sessiz") Mutasyon Etki Laboratuvarı',
    mutSub: "Aşağıdaki herhangi bir kodona tıklayıp aynı amino asidi kodlayan başka bir kodon seçin. Protein dizisi %100 aynı kalırken ribozom trafiğinin nasıl değiştiğini izleyin.",
    scanBtn: "En Tehlikeli Sessiz Mutasyonları Tara",
    resetMutBtn: "Değişiklikleri Sıfırla",
    hotspotsTitle: "En Yüksek Çarpışma Riski Yaratan Tek Nükleotidlik Sessiz Mutasyonlar",
    mutCompTitle: "Yabanıl Tip (WT) vs. Sessiz Mutant Karşılaştırması",
    proteinIdentical: "✓ Protein Dizisi %100 Birebir Aynı",
    loadMutantBtn: "Mutantı Canlı Simülatöre Yükle",
    harmTitle: "Duraklama Koruyucu Kodon Harmonizasyonu vs. Klasik CAI Optimizasyonu",
    harmSub: "Bir geni doğal organizmasından hedef konağa aktarırken katlanma duraklama ritmini (%MinMax korelasyonu) koruyun ve ribozom tıkanmalarını önleyin.",
    sourceOrgLabel: "Doğal Kaynak Organizma",
    targetOrgLabel: "Hedef Ekspresyon Konağı",
    runHarmBtn: "Karşılaştır ve Harmonize Et",
    thVariant: "Varyant",
    thRhythm: "Katlanma Ritmi (r)",
    thFlux: "Protein Akışı (/dk)",
    thColl: "Çarpışma %",
    thNgd: "NGD Riski",
    thChanged: "Değişen Kodon",
    harmChartTitle: "%MinMax Eş-Zamanlı Katlanma Ritmi Hizalaması",
    legNative: "Doğal Kaynaktaki Ritim (Hedef)",
    legHarm: "RiboFlow Harmonize Dizi",
    legNaive: "Klasik CAI Optimizasyonu (Düzleşmiş Duraklamalar)",
    harmFastaTitle: "Harmonize Edilmiş CDS Çıktısı (%100 Aynı Protein)",
    copyFastaBtn: "FASTA Kopyala",
    loadHarmBtn: "Harmonize Geni Simüle Et",
    pause: "Duraklat",
    play: "Oynat",
    copied: "Kopyalandı!",
  },
  en: {
    subtitle: "Co-Translational Ribosome Traffic (ℓ-TASEP) & Codon Harmonization Lab",
    sampleLabel: "Gene Preset / CDS",
    hostLabel: "Translation Host Organism",
    alphaLabel: "Initiation Rate (α)",
    runSim: "Run Simulation",
    customFasta: "Custom FASTA",
    pasteFasta: "Paste coding DNA / mRNA sequence (FASTA or raw ATG...TAA):",
    loadSeq: "Analyze & Simulate Sequence",
    kpiFlux: "Protein Flux (J)",
    kpiFluxSub: "proteins / min",
    kpiCollision: "Disome Collision Time",
    kpiNgd: "NGD / RQC Risk",
    kpiCai: "CAI & GC / GC3",
    kpiPhase: "TASEP Traffic Phase",
    tabTraffic: "1. Live Ribosome Traffic & Ribo-seq",
    tabMutation: "2. Silent Mutation Lab",
    tabHarmonize: "3. Codon Harmonizer",
    stageTitle: "Real-Time ℓ-TASEP Ribosome Translation Visualizer",
    stageHint: "Each ribosome occludes 10 codons (~30 nt). Red pulses mark steric disome collisions that trigger No-Go Decay (NGD).",
    speedLabel: "Speed",
    legFast: "Fast codon",
    legMid: "Moderate codon",
    legSlow: "Rare / pause codon",
    legRibo: "Translating ribosome",
    legColl: "Collided disome (NGD risk)",
    chartRiboseqTitle: "Simulated Ribo-seq & Disome-seq Profile",
    chartRiboseqSub: "Time-averaged ribosome A-site occupancy (teal) and disome collision peaks (crimson) along the mRNA.",
    chartMinmaxTitle: "%MinMax Co-Translational Folding Rhythm",
    chartMinmaxSub: "Positive peaks indicate fast clusters; negative amber troughs indicate putative co-translational folding pauses.",
    mutTitle: 'Synonymous ("Silent") Mutation Impact Lab',
    mutSub: "Click any codon below to swap it for a synonymous codon. The protein sequence stays 100% identical—watch how translation kinetics and collisions shift.",
    scanBtn: "Scan Most Dangerous Silent Mutations",
    resetMutBtn: "Reset Edits",
    hotspotsTitle: "Top Single-Nucleotide Silent Collision Triggers",
    mutCompTitle: "Wild-Type vs. Silent Mutant Comparison",
    proteinIdentical: "✓ Protein Sequence 100% Identical",
    loadMutantBtn: "Load Mutant into Live Canvas",
    harmTitle: "Pause-Preserving Codon Harmonizer vs. Naive CAI Optimization",
    harmSub: "Transfer a gene from its native organism to a recombinant host while preserving its co-translational folding rhythm (%MinMax correlation) and preventing disome bottlenecks.",
    sourceOrgLabel: "Native Source Organism",
    targetOrgLabel: "Target Expression Host",
    runHarmBtn: "Compare & Harmonize",
    thVariant: "Variant",
    thRhythm: "Folding Rhythm (r)",
    thFlux: "Protein Flux (/min)",
    thColl: "Collision %",
    thNgd: "NGD Risk",
    thChanged: "Codons Changed",
    harmChartTitle: "%MinMax Co-Translational Rhythm Alignment",
    legNative: "Native in Source (Target Rhythm)",
    legHarm: "RiboFlow Harmonized",
    legNaive: "Naive CAI Optimized (Flattened Pauses)",
    harmFastaTitle: "Harmonized CDS Output (100% Identical Protein)",
    copyFastaBtn: "Copy FASTA",
    loadHarmBtn: "Simulate Harmonized Gene",
    pause: "Pause",
    play: "Play",
    copied: "Copied!",
  },
};

const state = {
  lang: "tr",
  meta: null,
  currentSequence: "",
  currentName: "",
  organism: "ecoli",
  alpha: 0.85,
  analysis: null,
  simulation: null,
  playing: true,
  speed: 1.0,
  frameCursor: 0,
  lastTick: 0,
  edits: {},
  selectedCodon: null,
  lastMutComparison: null,
  lastHarmReport: null,
};

const $ = (id) => document.getElementById(id);

function applyLanguage() {
  const dict = I18N[state.lang];
  document.querySelectorAll("[data-i18n]").forEach((el) => {
    const key = el.getAttribute("data-i18n");
    if (dict[key]) el.textContent = dict[key];
  });
  $("play-pause-btn").textContent = state.playing ? dict.pause : dict.play;
  $("lang-tr").classList.toggle("active", state.lang === "tr");
  $("lang-en").classList.toggle("active", state.lang === "en");
  updateSampleDescription();
}

function updateSampleDescription() {
  if (!state.meta) return;
  const sampleId = $("sample-select").value;
  const s = state.meta.samples.find((x) => x.id === sampleId);
  if (s) {
    $("sample-description").textContent = state.lang === "tr" ? s.summary_tr : s.summary_en;
  }
}

async function apiPost(endpoint, payload) {
  const res = await fetch(endpoint, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(payload),
  });
  const data = await res.json();
  if (!res.ok) throw new Error(data.error || "API error");
  return data;
}

async function init() {
  const res = await fetch("/api/meta");
  state.meta = await res.json();

  // Populate samples
  const sampleSel = $("sample-select");
  state.meta.samples.forEach((s) => {
    const opt = document.createElement("option");
    opt.value = s.id;
    opt.textContent = s.name;
    sampleSel.appendChild(opt);
  });

  // Populate organism selects
  ["organism-select", "source-org-select", "target-org-select"].forEach((id) => {
    const sel = $(id);
    state.meta.organisms.forEach((o) => {
      const opt = document.createElement("option");
      opt.value = o.key;
      opt.textContent = `${o.name} (${o.common_name})`;
      sel.appendChild(opt);
    });
  });

  // Default to Ribosome Traffic-Jam Benchmark in E. coli for immediate visual punch,
  // or Human Preproinsulin in E. coli
  sampleSel.value = "human_ins";
  $("organism-select").value = "ecoli";
  $("source-org-select").value = "human";
  $("target-org-select").value = "ecoli";

  const firstSample = state.meta.samples.find((s) => s.id === "human_ins");
  state.currentSequence = firstSample.sequence;
  state.currentName = firstSample.name;

  bindEvents();
  applyLanguage();
  await runSimulation();
  requestAnimationFrame(animationLoop);
}

function bindEvents() {
  $("lang-tr").addEventListener("click", () => { state.lang = "tr"; applyLanguage(); });
  $("lang-en").addEventListener("click", () => { state.lang = "en"; applyLanguage(); });

  $("sample-select").addEventListener("change", async () => {
    const s = state.meta.samples.find((x) => x.id === $("sample-select").value);
    if (!s) return;
    state.currentSequence = s.sequence;
    state.currentName = s.name;
    state.edits = {};
    state.selectedCodon = null;
    $("organism-select").value = s.recommended_host;
    $("source-org-select").value = s.organism;
    $("target-org-select").value = s.recommended_host;
    const hostObj = state.meta.organisms.find((o) => o.key === s.recommended_host);
    if (hostObj) {
      state.alpha = hostObj.default_initiation_rate;
      $("alpha-slider").value = state.alpha;
      $("alpha-value").textContent = `${state.alpha.toFixed(2)} /s`;
    }
    updateSampleDescription();
    $("mutation-result-card").classList.add("hidden");
    $("hotspots-panel").classList.add("hidden");
    await runSimulation();
  });

  $("organism-select").addEventListener("change", async () => {
    state.organism = $("organism-select").value;
    const orgObj = state.meta.organisms.find((o) => o.key === state.organism);
    if (orgObj) {
      state.alpha = orgObj.default_initiation_rate;
      $("alpha-slider").value = state.alpha;
      $("alpha-value").textContent = `${state.alpha.toFixed(2)} /s`;
    }
    state.edits = {};
    $("mutation-result-card").classList.add("hidden");
    $("hotspots-panel").classList.add("hidden");
    await runSimulation();
  });

  $("alpha-slider").addEventListener("input", () => {
    state.alpha = parseFloat($("alpha-slider").value);
    $("alpha-value").textContent = `${state.alpha.toFixed(2)} /s`;
  });
  $("alpha-slider").addEventListener("change", () => runSimulation());

  $("run-btn").addEventListener("click", () => runSimulation());

  $("toggle-fasta-btn").addEventListener("click", () => {
    $("custom-fasta-box").classList.toggle("hidden");
    if (!$("fasta-input").value.trim()) {
      $("fasta-input").value = `>${state.currentName}\n${state.currentSequence}`;
    }
  });

  $("apply-fasta-btn").addEventListener("click", async () => {
    $("fasta-error").textContent = "";
    try {
      state.currentSequence = $("fasta-input").value;
      state.currentName = "Custom CDS";
      state.edits = {};
      await runSimulation();
      $("custom-fasta-box").classList.add("hidden");
    } catch (err) {
      $("fasta-error").textContent = err.message;
    }
  });

  // Tabs
  document.querySelectorAll(".tab-btn").forEach((btn) => {
    btn.addEventListener("click", () => {
      const target = btn.dataset.tab;
      document.querySelectorAll(".tab-btn").forEach((b) => b.classList.toggle("active", b === btn));
      document.querySelectorAll(".tab-panel").forEach((p) => p.classList.toggle("active", p.id === `tab-${target}`));
      if (target === "traffic") {
        drawRiboseqChart();
        drawMinmaxChart();
      }
    });
  });

  // Playback
  $("play-pause-btn").addEventListener("click", () => {
    state.playing = !state.playing;
    $("play-pause-btn").textContent = state.playing ? I18N[state.lang].pause : I18N[state.lang].play;
  });

  $("speed-select").addEventListener("change", () => {
    state.speed = parseFloat($("speed-select").value);
  });

  // Mutation lab buttons
  $("scan-hotspots-btn").addEventListener("click", scanHotspots);
  $("reset-mut-btn").addEventListener("click", () => {
    state.edits = {};
    state.selectedCodon = null;
    $("codon-editor-drawer").classList.add("hidden");
    $("mutation-result-card").classList.add("hidden");
    renderCodonStrip();
  });
  $("load-mutant-btn").addEventListener("click", async () => {
    if (!state.lastMutComparison) return;
    state.currentSequence = state.lastMutComparison.mutant_sequence;
    state.edits = {};
    $("mutation-result-card").classList.add("hidden");
    document.querySelector('[data-tab="traffic"]').click();
    await runSimulation();
  });

  // Harmonizer
  $("run-harm-btn").addEventListener("click", runHarmonization);
  $("copy-harm-btn").addEventListener("click", async () => {
    if (!state.lastHarmReport) return;
    await navigator.clipboard.writeText($("harm-fasta-out").textContent);
    const btn = $("copy-harm-btn");
    const prev = btn.textContent;
    btn.textContent = I18N[state.lang].copied;
    setTimeout(() => { btn.textContent = prev; }, 1400);
  });
  $("load-harm-btn").addEventListener("click", async () => {
    if (!state.lastHarmReport) return;
    state.currentSequence = state.lastHarmReport.harmonized.sequence;
    $("organism-select").value = state.lastHarmReport.host_organism;
    state.organism = state.lastHarmReport.host_organism;
    state.edits = {};
    document.querySelector('[data-tab="traffic"]').click();
    await runSimulation();
  });
}

async function runSimulation() {
  state.organism = $("organism-select").value;
  state.alpha = parseFloat($("alpha-slider").value);

  const data = await apiPost("/api/simulate", {
    sequence: state.currentSequence,
    organism: state.organism,
    initiation_rate: state.alpha,
    duration: 120,
    seed: 42,
  });

  state.analysis = data.analysis;
  state.simulation = data.simulation;
  state.currentSequence = data.analysis.sequence;
  state.frameCursor = 0;

  updateKPIs();
  drawRiboseqChart();
  drawMinmaxChart();
  renderCodonStrip();
}

function updateKPIs() {
  const sim = state.simulation;
  const an = state.analysis;

  $("kpi-flux").textContent = sim.protein_flux_per_min.toFixed(1);
  $("kpi-collision").textContent = `${(sim.collision_fraction * 100).toFixed(1)}%`;
  $("kpi-collision-events").textContent = `${sim.collision_events} events • bottleneck #${sim.bottleneck_codon + 1}`;

  const ngdEl = $("kpi-ngd");
  ngdEl.textContent = `${sim.ngd_risk_score.toFixed(1)} / 100`;
  ngdEl.className = "kpi-val " + (sim.ngd_risk_score >= 45 ? "risk-high" : sim.ngd_risk_score >= 18 ? "risk-mid" : "risk-low");
  $("kpi-queue").textContent = `max queue: ${sim.max_queue_length} ribos`;

  $("kpi-cai").textContent = an.cai.toFixed(3);
  $("kpi-gc").textContent = `GC ${an.gc_percent}% • GC3 ${an.gc3_percent}%`;

  $("kpi-phase").textContent = sim.phase_regime;
  $("kpi-ribos").textContent = `avg ${sim.mean_ribosomes_on_mrna} ribos on ${an.codon_count} codons`;
}

function codonSpeedColor(rate, maxRate) {
  const norm = Math.min(1, Math.max(0, rate / Math.max(maxRate, 1)));
  if (norm < 0.35) return "#f43f5e";
  if (norm < 0.62) return "#f59e0b";
  return "#10b981";
}

function animationLoop(ts) {
  if (!state.lastTick) state.lastTick = ts;
  const dt = (ts - state.lastTick) / 1000;
  state.lastTick = ts;

  if (state.simulation && state.simulation.frames.length > 1) {
    if (state.playing) {
      state.frameCursor += dt * 7.5 * state.speed;
      if (state.frameCursor >= state.simulation.frames.length - 1) {
        state.frameCursor = 0;
      }
    }
    drawTrafficCanvas();
  }
  requestAnimationFrame(animationLoop);
}

function drawTrafficCanvas() {
  const canvas = $("traffic-canvas");
  const ctx = canvas.getContext("2d");
  const W = canvas.width;
  const H = canvas.height;
  const sim = state.simulation;
  const frames = sim.frames;
  if (!frames || !frames.length) return;

  const idxA = Math.floor(state.frameCursor);
  const idxB = Math.min(frames.length - 1, idxA + 1);
  const frac = state.frameCursor - idxA;
  const frameA = frames[idxA];
  const frameB = frames[idxB];

  $("sim-clock").textContent = `t = ${frameA.time.toFixed(1)} s • ${frameA.completed} proteins`;

  ctx.clearRect(0, 0, W, H);

  // Subtle background grid
  ctx.strokeStyle = "rgba(35, 50, 77, 0.35)";
  ctx.lineWidth = 1;
  for (let x = 50; x < W - 40; x += 80) {
    ctx.beginPath();
    ctx.moveTo(x, 20);
    ctx.lineTo(x, H - 25);
    ctx.stroke();
  }

  const padL = 54;
  const padR = 54;
  const trackW = W - padL - padR;
  const trackY = 178;
  const L = sim.codon_count;
  const codonW = trackW / L;
  const maxRate = Math.max(...sim.codon_rates);

  // 5' and 3' labels
  ctx.fillStyle = "#8ea3c4";
  ctx.font = "bold 12px monospace";
  ctx.fillText("5' AUG", 10, trackY + 4);
  ctx.fillText("3' STOP", W - 49, trackY + 4);

  // Draw mRNA codon speed bar
  for (let i = 0; i < L; i++) {
    const x = padL + i * codonW;
    ctx.fillStyle = codonSpeedColor(sim.codon_rates[i], maxRate);
    ctx.fillRect(x, trackY - 5, Math.max(1.5, codonW - 0.6), 10);
  }

  // Bottleneck marker
  const bx = padL + (sim.bottleneck_codon + 0.5) * codonW;
  ctx.strokeStyle = "rgba(244, 63, 94, 0.55)";
  ctx.setLineDash([4, 4]);
  ctx.beginPath();
  ctx.moveTo(bx, 30);
  ctx.lineTo(bx, trackY + 38);
  ctx.stroke();
  ctx.setLineDash([]);
  ctx.fillStyle = "#fda4af";
  ctx.font = "11px monospace";
  ctx.fillText(`Bottleneck #${sim.bottleneck_codon + 1}`, Math.min(W - 140, Math.max(padL, bx - 42)), trackY + 54);

  // Interpolate ribosomes by ID between frameA and frameB
  const nextMap = new Map(frameB.ribosomes.map((r) => [r.id, r]));

  frameA.ribosomes.forEach((r) => {
    const next = nextMap.get(r.id);
    const interpPos = next ? r.pos + (next.pos - r.pos) * frac : r.pos;
    const cx = padL + (interpPos + 0.5) * codonW;
    const colliding = r.colliding;

    // Footprint span (10 codons)
    const fpW = Math.max(18, sim.ribosome_footprint * codonW);
    ctx.fillStyle = colliding ? "rgba(244, 63, 94, 0.22)" : "rgba(45, 212, 191, 0.12)";
    ctx.fillRect(cx - fpW * 0.5, trackY - 26, fpW, 52);

    // Growing nascent polypeptide chain (beads rising above ribosome)
    const progress = interpPos / Math.max(1, L - 1);
    const beads = Math.max(2, Math.floor(progress * 11));
    ctx.strokeStyle = colliding ? "#fda4af" : "#7dd3fc";
    ctx.lineWidth = 2;
    ctx.beginPath();
    let px = cx;
    let py = trackY - 26;
    ctx.moveTo(px, py);
    for (let b = 1; b <= beads; b++) {
      px += Math.sin(b * 1.3 + r.id) * 6;
      py -= 9;
      ctx.lineTo(px, py);
    }
    ctx.stroke();

    // Ribosome large subunit (top) & small subunit (bottom)
    ctx.fillStyle = colliding ? "#f43f5e" : "#2dd4bf";
    ctx.beginPath();
    ctx.ellipse(cx, trackY - 15, Math.min(20, fpW * 0.44), 12, 0, 0, Math.PI * 2);
    ctx.fill();

    ctx.fillStyle = colliding ? "#9f1239" : "#0f766e";
    ctx.beginPath();
    ctx.ellipse(cx, trackY + 11, Math.min(15, fpW * 0.34), 7, 0, 0, Math.PI * 2);
    ctx.fill();

    // Disome collision pulse ring
    if (colliding) {
      const pulse = 18 + Math.sin(performance.now() / 120) * 5;
      ctx.strokeStyle = "rgba(255, 46, 99, 0.85)";
      ctx.lineWidth = 2.2;
      ctx.beginPath();
      ctx.arc(cx + fpW * 0.45, trackY - 4, pulse, 0, Math.PI * 2);
      ctx.stroke();
    }
  });
}

function drawRiboseqChart() {
  const canvas = $("riboseq-canvas");
  const ctx = canvas.getContext("2d");
  const sim = state.simulation;
  if (!sim) return;

  const W = canvas.width;
  const H = canvas.height;
  ctx.clearRect(0, 0, W, H);

  const pad = { l: 42, r: 16, t: 18, b: 28 };
  const plotW = W - pad.l - pad.r;
  const plotH = H - pad.t - pad.b;
  const L = sim.codon_occupancy.length;
  const maxY = Math.max(0.12, ...sim.codon_occupancy, ...sim.collision_profile);

  // Axes
  ctx.strokeStyle = "#23324d";
  ctx.beginPath();
  ctx.moveTo(pad.l, pad.t);
  ctx.lineTo(pad.l, H - pad.b);
  ctx.lineTo(W - pad.r, H - pad.b);
  ctx.stroke();

  // Ribo-seq area
  ctx.fillStyle = "rgba(45, 212, 191, 0.25)";
  ctx.strokeStyle = "#2dd4bf";
  ctx.lineWidth = 1.8;
  ctx.beginPath();
  ctx.moveTo(pad.l, H - pad.b);
  for (let i = 0; i < L; i++) {
    const x = pad.l + (i / Math.max(1, L - 1)) * plotW;
    const y = H - pad.b - (sim.codon_occupancy[i] / maxY) * plotH;
    ctx.lineTo(x, y);
  }
  ctx.lineTo(W - pad.r, H - pad.b);
  ctx.closePath();
  ctx.fill();
  ctx.stroke();

  // Disome-seq overlay
  ctx.fillStyle = "rgba(244, 63, 94, 0.45)";
  ctx.strokeStyle = "#f43f5e";
  ctx.lineWidth = 1.8;
  ctx.beginPath();
  ctx.moveTo(pad.l, H - pad.b);
  for (let i = 0; i < L; i++) {
    const x = pad.l + (i / Math.max(1, L - 1)) * plotW;
    const y = H - pad.b - (sim.collision_profile[i] / maxY) * plotH;
    ctx.lineTo(x, y);
  }
  ctx.lineTo(W - pad.r, H - pad.b);
  ctx.closePath();
  ctx.fill();
  ctx.stroke();

  ctx.fillStyle = "#8ea3c4";
  ctx.font = "11px monospace";
  ctx.fillText(maxY.toFixed(2), 4, pad.t + 8);
  ctx.fillText("0.00", 8, H - pad.b);
  ctx.fillText("Codon 1", pad.l, H - 8);
  ctx.fillText(`Codon ${L}`, W - 74, H - 8);
}

function drawMinmaxChart() {
  const canvas = $("minmax-canvas");
  const ctx = canvas.getContext("2d");
  const an = state.analysis;
  if (!an) return;

  const W = canvas.width;
  const H = canvas.height;
  ctx.clearRect(0, 0, W, H);

  const pad = { l: 42, r: 16, t: 18, b: 28 };
  const plotW = W - pad.l - pad.r;
  const plotH = H - pad.t - pad.b;
  const midY = pad.t + plotH / 2;
  const scores = an.codons.map((c) => c.minmax_score);
  const L = scores.length;

  // Zero axis
  ctx.strokeStyle = "#334155";
  ctx.setLineDash([4, 4]);
  ctx.beginPath();
  ctx.moveTo(pad.l, midY);
  ctx.lineTo(W - pad.r, midY);
  ctx.stroke();
  ctx.setLineDash([]);

  // Curve
  for (let i = 0; i < L - 1; i++) {
    const x1 = pad.l + (i / (L - 1)) * plotW;
    const y1 = midY - (scores[i] / 100) * (plotH / 2);
    const x2 = pad.l + ((i + 1) / (L - 1)) * plotW;
    const y2 = midY - (scores[i + 1] / 100) * (plotH / 2);
    ctx.strokeStyle = scores[i] < -10 ? "#f59e0b" : "#38bdf8";
    ctx.lineWidth = 2.2;
    ctx.beginPath();
    ctx.moveTo(x1, y1);
    ctx.lineTo(x2, y2);
    ctx.stroke();
  }

  ctx.fillStyle = "#8ea3c4";
  ctx.font = "11px monospace";
  ctx.fillText("+100", 4, pad.t + 8);
  ctx.fillText("   0", 4, midY + 4);
  ctx.fillText("-100", 4, H - pad.b);
  ctx.fillText("Fast (+)", pad.l + 6, pad.t + 12);
  ctx.fillText("Pause (-)", pad.l + 6, H - pad.b - 6);
}

/* ---------- TAB 2: SILENT MUTATION LAB ---------- */

function renderCodonStrip() {
  const strip = $("codon-strip");
  strip.innerHTML = "";
  if (!state.analysis) return;

  const codons = state.analysis.codons;
  const maxRate = Math.max(...codons.map((c) => c.elongation_rate));

  codons.forEach((c) => {
    const locked = c.index === 0 || c.amino_acid === "*" || c.synonymous_options.length <= 1;
    const activeCodon = state.edits[c.index] || c.codon;
    const btn = document.createElement("button");
    btn.className = "codon-pill";
    if (locked) btn.classList.add("locked");
    if (state.edits[c.index] && state.edits[c.index] !== c.codon) btn.classList.add("edited");
    if (state.selectedCodon === c.index) btn.classList.add("selected");

    btn.style.borderBottom = `3px solid ${codonSpeedColor(c.elongation_rate, maxRate)}`;
    btn.innerHTML = `<span class="codon-idx">${c.index + 1}</span><span class="codon-aa">${c.amino_acid}</span><span>${activeCodon}</span>`;

    if (!locked) {
      btn.addEventListener("click", () => openCodonEditor(c));
    }
    strip.appendChild(btn);
  });
}

function openCodonEditor(codonObj) {
  state.selectedCodon = codonObj.index;
  renderCodonStrip();

  const drawer = $("codon-editor-drawer");
  drawer.classList.remove("hidden");
  $("editor-title").textContent = `Codon #${codonObj.index + 1} — ${codonObj.aa_three} (${codonObj.amino_acid})`;
  $("editor-meta").textContent = `Wild-Type: ${codonObj.codon} (${codonObj.elongation_rate} codons/s)`;

  const container = $("synonymous-buttons");
  container.innerHTML = "";
  const activeCodon = state.edits[codonObj.index] || codonObj.codon;

  codonObj.synonymous_options.forEach((syn) => {
    const btn = document.createElement("button");
    btn.className = "syn-btn" + (syn === activeCodon ? " current" : "");
    btn.textContent = syn === codonObj.codon ? `${syn} (WT)` : syn;
    btn.addEventListener("click", async () => {
      if (syn === codonObj.codon) {
        delete state.edits[codonObj.index];
      } else {
        state.edits[codonObj.index] = syn;
      }
      renderCodonStrip();
      openCodonEditor(codonObj);
      if (Object.keys(state.edits).length > 0) {
        await runMutationComparison();
      } else {
        $("mutation-result-card").classList.add("hidden");
      }
    });
    container.appendChild(btn);
  });
}

async function runMutationComparison() {
  const comp = await apiPost("/api/mutate", {
    sequence: state.currentSequence,
    organism: state.organism,
    edits: state.edits,
    initiation_rate: state.alpha,
    seed: 42,
  });
  state.lastMutComparison = comp;

  $("mutation-result-card").classList.remove("hidden");
  const metricsEl = $("mut-metrics");
  const fluxSign = comp.delta_flux_per_min >= 0 ? "+" : "";
  const ngdSign = comp.delta_ngd_risk >= 0 ? "+" : "";
  const collDeltaPct = (comp.delta_collision_fraction * 100).toFixed(2);

  metricsEl.innerHTML = `
    <div class="mut-stat">
      <span>Protein Flux (J)</span>
      <strong>${comp.wildtype_sim.protein_flux_per_min} → ${comp.mutant_sim.protein_flux_per_min} (${fluxSign}${comp.delta_flux_percent}%)</strong>
    </div>
    <div class="mut-stat">
      <span>Disome Collision %</span>
      <strong>${(comp.wildtype_sim.collision_fraction * 100).toFixed(1)}% → ${(comp.mutant_sim.collision_fraction * 100).toFixed(1)}% (${collDeltaPct >= 0 ? "+" : ""}${collDeltaPct}%)</strong>
    </div>
    <div class="mut-stat">
      <span>NGD / RQC Risk</span>
      <strong class="${comp.delta_ngd_risk > 5 ? "risk-high" : "risk-low"}">${comp.wildtype_sim.ngd_risk_score} → ${comp.mutant_sim.ngd_risk_score} (${ngdSign}${comp.delta_ngd_risk})</strong>
    </div>
    <div class="mut-stat">
      <span>CAI Change</span>
      <strong>${comp.wildtype_cai.toFixed(3)} → ${comp.mutant_cai.toFixed(3)}</strong>
    </div>
  `;

  drawMutationDiffCanvas(comp);
}

function drawMutationDiffCanvas(comp) {
  const canvas = $("mut-diff-canvas");
  const ctx = canvas.getContext("2d");
  const W = canvas.width;
  const H = canvas.height;
  ctx.clearRect(0, 0, W, H);

  const pad = { l: 48, r: 20, t: 22, b: 26 };
  const plotW = W - pad.l - pad.r;
  const plotH = H - pad.t - pad.b;

  const wt = comp.wildtype_sim.codon_occupancy;
  const mut = comp.mutant_sim.codon_occupancy;
  const L = wt.length;
  const maxY = Math.max(0.12, ...wt, ...mut);

  ctx.strokeStyle = "#23324d";
  ctx.beginPath();
  ctx.moveTo(pad.l, pad.t);
  ctx.lineTo(pad.l, H - pad.b);
  ctx.lineTo(W - pad.r, H - pad.b);
  ctx.stroke();

  const drawLine = (arr, color, dashed) => {
    ctx.strokeStyle = color;
    ctx.lineWidth = 2;
    ctx.setLineDash(dashed ? [5, 4] : []);
    ctx.beginPath();
    arr.forEach((v, i) => {
      const x = pad.l + (i / (L - 1)) * plotW;
      const y = H - pad.b - (v / maxY) * plotH;
      if (i === 0) ctx.moveTo(x, y);
      else ctx.lineTo(x, y);
    });
    ctx.stroke();
    ctx.setLineDash([]);
  };

  drawLine(wt, "#94a3b8", true);
  drawLine(mut, "#f43f5e", false);

  // Mark edited codons
  comp.edits.forEach((e) => {
    const x = pad.l + (e.index / (L - 1)) * plotW;
    ctx.strokeStyle = "#f59e0b";
    ctx.beginPath();
    ctx.moveTo(x, pad.t);
    ctx.lineTo(x, H - pad.b);
    ctx.stroke();
    ctx.fillStyle = "#fbbf24";
    ctx.font = "11px monospace";
    ctx.fillText(`#${e.index + 1} ${e.from_codon}→${e.to_codon}`, Math.min(W - 140, x + 4), pad.t + 12);
  });

  ctx.fillStyle = "#94a3b8";
  ctx.font = "11px monospace";
  ctx.fillText("Dashed: Wild-Type Ribo-seq | Red: Silent Mutant Ribo-seq", pad.l + 6, H - 6);
}

async function scanHotspots() {
  const btn = $("scan-hotspots-btn");
  btn.disabled = true;
  try {
    const data = await apiPost("/api/scan-mutations", {
      sequence: state.currentSequence,
      organism: state.organism,
      initiation_rate: state.alpha,
      top_k: 6,
      seed: 42,
    });
    const panel = $("hotspots-panel");
    const list = $("hotspots-list");
    list.innerHTML = "";
    panel.classList.remove("hidden");

    data.hotspots.forEach((h) => {
      const el = document.createElement("div");
      el.className = "hotspot-item";
      el.innerHTML = `
        <strong>#${h.index + 1} ${h.aa_three}: ${h.from_codon} → ${h.to_codon}</strong>
        <div class="muted">Rate -${h.rate_drop_percent}% • ΔNGD <span class="risk-high">${h.delta_ngd_risk >= 0 ? "+" : ""}${h.delta_ngd_risk}</span> • ΔFlux ${h.delta_flux_per_min >= 0 ? "+" : ""}${h.delta_flux_per_min}</div>
      `;
      el.addEventListener("click", async () => {
        state.edits = { [h.index]: h.to_codon };
        renderCodonStrip();
        await runMutationComparison();
      });
      list.appendChild(el);
    });
  } finally {
    btn.disabled = false;
  }
}

/* ---------- TAB 3: CODON HARMONIZER ---------- */

async function runHarmonization() {
  const btn = $("run-harm-btn");
  btn.disabled = true;
  try {
    const rep = await apiPost("/api/harmonize", {
      sequence: state.currentSequence,
      source: $("source-org-select").value,
      host: $("target-org-select").value,
      initiation_rate: state.alpha,
      seed: 42,
    });
    state.lastHarmReport = rep;
    $("harm-results").classList.remove("hidden");

    const tbody = $("harm-tbody");
    tbody.innerHTML = "";
    const variants = [rep.wildtype_in_host, rep.naive_optimized, rep.harmonized];
    variants.forEach((v, idx) => {
      const tr = document.createElement("tr");
      if (idx === 2) tr.className = "highlight-row";
      tr.innerHTML = `
        <td><strong>${v.label}</strong></td>
        <td class="mono-cell">${v.cai_in_host.toFixed(3)}</td>
        <td class="mono-cell">${v.rhythm_correlation >= 0 ? "+" : ""}${v.rhythm_correlation.toFixed(3)}</td>
        <td class="mono-cell">${v.simulation.protein_flux_per_min.toFixed(2)}</td>
        <td class="mono-cell">${(v.simulation.collision_fraction * 100).toFixed(2)}%</td>
        <td class="mono-cell">${v.simulation.ngd_risk_score.toFixed(1)}</td>
        <td class="mono-cell">${v.codons_changed} (${v.percent_changed}%)</td>
      `;
      tbody.appendChild(tr);
    });

    drawHarmonizationCanvas(rep);
    $("harm-fasta-out").textContent =
      `>${state.currentName}_harmonized_${rep.source_organism}_to_${rep.host_organism}\n${rep.harmonized.sequence}`;
  } finally {
    btn.disabled = false;
  }
}

function drawHarmonizationCanvas(rep) {
  const canvas = $("harm-canvas");
  const ctx = canvas.getContext("2d");
  const W = canvas.width;
  const H = canvas.height;
  ctx.clearRect(0, 0, W, H);

  const pad = { l: 46, r: 18, t: 18, b: 26 };
  const plotW = W - pad.l - pad.r;
  const plotH = H - pad.t - pad.b;
  const midY = pad.t + plotH / 2;

  ctx.strokeStyle = "#334155";
  ctx.setLineDash([4, 4]);
  ctx.beginPath();
  ctx.moveTo(pad.l, midY);
  ctx.lineTo(W - pad.r, midY);
  ctx.stroke();
  ctx.setLineDash([]);

  const drawSeries = (arr, color, width, dashed) => {
    ctx.strokeStyle = color;
    ctx.lineWidth = width;
    ctx.setLineDash(dashed ? [6, 4] : []);
    ctx.beginPath();
    arr.forEach((v, i) => {
      const x = pad.l + (i / Math.max(1, arr.length - 1)) * plotW;
      const y = midY - (v / 100) * (plotH / 2);
      if (i === 0) ctx.moveTo(x, y);
      else ctx.lineTo(x, y);
    });
    ctx.stroke();
    ctx.setLineDash([]);
  };

  drawSeries(rep.naive_optimized.minmax_profile, "#a78bfa", 1.8, false);
  drawSeries(rep.native_source_minmax, "#e2e8f0", 2.2, true);
  drawSeries(rep.harmonized.minmax_profile, "#2dd4bf", 2.4, false);

  ctx.fillStyle = "#8ea3c4";
  ctx.font = "11px monospace";
  ctx.fillText("+100", 6, pad.t + 8);
  ctx.fillText("   0", 6, midY + 4);
  ctx.fillText("-100", 6, H - pad.b);
}

init();
