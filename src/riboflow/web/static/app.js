/* RiboFlow v2.0 Interactive Biophysical & Synthetic Biology Web Laboratory */

const I18N = {
  tr: {
    subtitle: "Biyofiziksel Ribozom Trafiği (ℓ-TASEP), 5' mRNA Katlanması & mRNA Aşı / SynBio Laboratuvarı",
    exportReport: "📄 Lab Raporu",
    sampleLabel: "Hazır Biyolojik Gen / CDS",
    hostLabel: "Çeviri Konak Organizması",
    alphaLabel: "Nominal Başlatma Hızı (α)",
    runSim: "Simüle Et",
    ncbiBtn: "🌐 NCBI Çek",
    customFasta: "FASTA",
    toggle5p: "5' mRNA Saç Tokası ΔG → α_eff Bağlantısı",
    toggleTunnel: "Ribozom Çıkış Tüneli Yük Freni (+K/R)",
    ncbiPrompt: "NCBI RefSeq'ten canlı kodlayan dizi (CDS) çek (örn: INS, TP53, CFTR, BRCA1, NM_000207):",
    ncbiFetchBtn: "Canlı CDS Çek",
    pasteFasta: "Kodlayan DNA / mRNA dizisini yapıştırın (FASTA veya ham ATG...TAA):",
    loadSeq: "Diziyi Analiz Et ve Simüle Et",
    kpiFlux: "Protein Akışı (J)",
    kpiHairpin: "5' Saç Tokası (ΔG)",
    kpiCollision: "Disome Çarpışma",
    kpiNgd: "NGD / RQC Riski",
    kpiTunnel: "Çıkış Tüneli Yükü",
    kpiPhase: "TASEP Fazı",
    tabTraffic: "1. Canlı Trafik & Tünel",
    tabFolding: "2. 5' mRNA Saç Tokası ΔG",
    tabMutation: "3. Sessiz Mutasyon Labı",
    tabSynbio: "4. SynBio & mRNA Aşı",
    tabPhase: "5. Faz Diyagramı (α Tarama)",
    stageTitle: "Gerçek Zamanlı ℓ-TASEP Ribozom Trafiği & Çıkış Tüneli Görselleştiricisi",
    stageHint: "Her ribozom 10 kodon (~30 nt) kaplar. Kırmızı halkalar = Disome çarpışması (NGD). Mor parlama = polibazik (+K/R) çıkış tüneli frenlemesi.",
    speedLabel: "Hız",
    legFast: "Hızlı kodon",
    legMid: "Orta kodon",
    legSlow: "Nadir / duraklama",
    legRibo: "Aktif ribozom",
    legTunnel: "Polibazik tünel freni (+K/R)",
    legColl: "Çarpışan disome (NGD)",
    chartRiboseqTitle: "Simüle Ribo-seq & Disome-seq",
    chartRiboseqSub: "Ribozom yoğunluğu (turkuaz) & çarpışma tepeleri (kırmızı).",
    chartMinmaxTitle: "%MinMax Katlanma Ritmi",
    chartMinmaxSub: "Pozitif hızlı bölgeler vs. negatif amber katlanma duraklamaları.",
    chartTunnelTitle: "Çıkış Tüneli Yükü (+K/R Freni)",
    chartTunnelSub: "Negatif rRNA tünelinde pozitif yüklü peptid elektrostatik frenlemesi.",
    foldTitle: "5' mRNA İkincil Yapısı (Turner MFE ΔG) & Saç Tokası Açıcı",
    foldSub: "AUG başlangıç kodonu etrafındaki ilk 45 nt'de oluşan sıkı saç tokaları ribozom bağlanmasını fiziksel olarak bloke eder. Tek tıkla proteini değiştirmeden saç tokasını eritin!",
    meltBtn: "🔥 5' Saç Tokasını Sinonim Olarak Erit",
    meltAfterTitle: "Sinonim Saç Tokası Eritme Sonrası (%100 Aynı Protein)",
    loadMeltedBtn: "Açılmış mRNA'yı Simülatöre Yükle",
    mutTitle: 'Sinonim ("Sessiz") Mutasyon Etki Laboratuvarı',
    mutSub: "Herhangi bir kodona tıklayıp aynı amino asidi kodlayan başka bir kodon seçin. Protein %100 aynı kalırken trafiğin nasıl değiştiğini izleyin.",
    scanBtn: "En Tehlikeli Sessiz Mutasyonları Tara",
    resetMutBtn: "Sıfırla",
    hotspotsTitle: "En Yüksek Çarpışma Riski Yaratan Tek Nükleotidlik Sessiz Mutasyonlar",
    mutCompTitle: "Yabanıl Tip (WT) vs. Sessiz Mutant Karşılaştırması",
    proteinIdentical: "✓ Protein Dizisi %100 Birebir Aynı",
    loadMutantBtn: "Mutantı Canlı Simülatöre Yükle",
    synTitle: "Çok Amaçlı Sentetik Biyoloji & mRNA Aşı Tasarımcısı",
    synSub: "Katlanma ritmini koruma, 5' saç tokasını açma, mRNA aşısı Üridin azaltma (TLR7/8 kaçışı), klonlama restriksiyon bölgelerini (EcoRI/BamHI/BsaI) temizleme ve homopolimer gidermeyi tek seferde optimize edin!",
    sourceOrgLabel: "Doğal Kaynak Organizma",
    targetOrgLabel: "Hedef Ekspresyon Konağı",
    runSynBtn: "🚀 SynBio Dizisini Tasarla",
    optPauses: "Katlanma Duraklamalarını Koru (%MinMax)",
    optUnfold: "5' Saç Tokasını Erit (ΔG_5')",
    optVaccine: "💉 mRNA Aşı Modu: Üridin Azaltma (m¹Ψ)",
    optScrub: "Klonlama Kesim Bölgelerini Temizle (EcoRI/BamHI/BsaI…)",
    optHp: "Homopolimerleri Temizle (≥5 nt)",
    thVariant: "Tasarım Stratejisi",
    thRhythm: "Katlanma Ritmi (r)",
    thFlux: "Akış (/dk)",
    thColl: "Çarpışma %",
    thNgd: "NGD Riski",
    harmChartTitle: "%MinMax Eş-Zamanlı Katlanma Ritmi Hizalaması",
    legNative: "Doğal Kaynaktaki Ritim (Hedef)",
    legSyn: "RiboFlow v2.0 SynBio Tasarımı",
    legNaive: "Klasik CAI Optimizasyonu (Düzleşmiş)",
    harmFastaTitle: "Senteze Hazır SynBio CDS FASTA (%100 Aynı Protein)",
    copyFastaBtn: "FASTA Kopyala",
    loadHarmBtn: "Tasarlanan Geni Simülatöre Yükle",
    phaseTitle: "TASEP İstatistiksel Fizik Faz Geçişi Diyagramı (α Taraması)",
    phaseSub: "Başlatma hızını α = 0.08 → 2.20 /s aralığında tarayarak mRNA'nın Düşük Yoğunluktan Trafik Tıkanıklığına geçtiği kritik eşik noktasını (α_c) hesaplar.",
    runPhaseBtn: "📈 Faz Geçişi Taramasını Başlat",
    legFluxWT: "Protein Akışı J(α) — Mevcut Gen",
    legNgdWT: "NGD Çarpışma Riski — Mevcut Gen",
    legFluxOpt: "Protein Akışı J(α) — SynBio Tasarımı",
    pause: "Duraklat",
    play: "Oynat",
    copied: "Kopyalandı!",
  },
  en: {
    subtitle: "Biophysical Ribosome Traffic (ℓ-TASEP), 5' mRNA Folding & mRNA Vaccine / SynBio Lab",
    exportReport: "📄 Lab Report",
    sampleLabel: "Gene Preset / CDS",
    hostLabel: "Translation Host Organism",
    alphaLabel: "Nominal Initiation Rate (α)",
    runSim: "Run Simulation",
    ncbiBtn: "🌐 NCBI",
    customFasta: "FASTA",
    toggle5p: "Couple 5' mRNA Hairpin ΔG → α_eff",
    toggleTunnel: "Ribosome Exit-Tunnel Charge Drag (+K/R)",
    ncbiPrompt: "Fetch live coding sequence from NCBI RefSeq (e.g. INS, TP53, CFTR, BRCA1, NM_000207):",
    ncbiFetchBtn: "Fetch Live CDS",
    pasteFasta: "Paste coding DNA / mRNA sequence (FASTA or raw ATG...TAA):",
    loadSeq: "Analyze & Simulate Sequence",
    kpiFlux: "Protein Flux (J)",
    kpiHairpin: "5' Hairpin (ΔG)",
    kpiCollision: "Disome Collisions",
    kpiNgd: "NGD / RQC Risk",
    kpiTunnel: "Exit-Tunnel Charge",
    kpiPhase: "TASEP Phase",
    tabTraffic: "1. Live Traffic & Tunnel",
    tabFolding: "2. 5' mRNA Hairpin ΔG",
    tabMutation: "3. Silent Mutation Lab",
    tabSynbio: "4. SynBio & mRNA Vaccine",
    tabPhase: "5. Phase Diagram (α Sweep)",
    stageTitle: "Real-Time ℓ-TASEP Ribosome Translation & Exit-Tunnel Visualizer",
    stageHint: "Each ribosome occludes 10 codons (~30 nt). Crimson pulses = Disome collisions (NGD). Violet glow = polybasic (+K/R) exit-tunnel drag.",
    speedLabel: "Speed",
    legFast: "Fast codon",
    legMid: "Moderate codon",
    legSlow: "Rare / pause",
    legRibo: "Active ribosome",
    legTunnel: "Polybasic exit-tunnel drag (+K/R)",
    legColl: "Collided disome (NGD)",
    chartRiboseqTitle: "Simulated Ribo-seq & Disome-seq",
    chartRiboseqSub: "Ribosome occupancy (teal) & disome collision peaks (crimson).",
    chartMinmaxTitle: "%MinMax Folding Rhythm",
    chartMinmaxSub: "Positive fast clusters vs. negative amber co-translational pauses.",
    chartTunnelTitle: "Exit-Tunnel Charge (+K/R Drag)",
    chartTunnelSub: "Positively charged nascent chain electrostatic braking inside rRNA tunnel.",
    foldTitle: "5' mRNA Secondary Structure (Turner MFE ΔG) & Hairpin Melter",
    foldSub: "Tight hairpins within the first 45 nt around AUG physically occlude 30S/43S recruitment. Click below to melt inhibitory hairpins synonymously (100% identical protein).",
    meltBtn: "🔥 Melt 5' Hairpin Synonymously",
    meltAfterTitle: "After Synonymous 5' Hairpin Melting (100% Identical Protein)",
    loadMeltedBtn: "Load Unfolded mRNA into Simulator",
    mutTitle: 'Synonymous ("Silent") Mutation Impact Lab',
    mutSub: "Click any codon below to swap it for a synonymous codon. Protein stays 100% identical—watch kinetics and collisions shift.",
    scanBtn: "Scan Most Dangerous Silent Mutations",
    resetMutBtn: "Reset",
    hotspotsTitle: "Top Single-Nucleotide Silent Collision Triggers",
    mutCompTitle: "Wild-Type vs. Silent Mutant Comparison",
    proteinIdentical: "✓ Protein Sequence 100% Identical",
    loadMutantBtn: "Load Mutant into Live Canvas",
    synTitle: "Multi-Objective SynBio & mRNA Vaccine Sequence Designer",
    synSub: "Co-optimize folding pause rhythm, 5' hairpin unfolding, mRNA vaccine Uridine depletion (TLR7/8 evasion), Golden Gate / cloning restriction-site scrubbing, and homopolymer removal—with 100% identical protein!",
    sourceOrgLabel: "Native Source Organism",
    targetOrgLabel: "Target Expression Host",
    runSynBtn: "🚀 Optimize SynBio Sequence",
    optPauses: "Preserve Folding Pauses (%MinMax)",
    optUnfold: "Melt 5' Hairpin (ΔG_5')",
    optVaccine: "💉 mRNA Vaccine Mode: Uridine Depletion (m¹Ψ)",
    optScrub: "Scrub Cloning Sites (EcoRI/BamHI/BsaI…)",
    optHp: "Remove Homopolymers (≥5 nt)",
    thVariant: "Design Strategy",
    thRhythm: "Folding Rhythm (r)",
    thFlux: "Flux (/min)",
    thColl: "Collision %",
    thNgd: "NGD Risk",
    harmChartTitle: "%MinMax Co-Translational Rhythm Alignment",
    legNative: "Native in Source (Target Rhythm)",
    legSyn: "RiboFlow v2.0 SynBio Designed",
    legNaive: "Naive CAI Optimized (Flattened)",
    harmFastaTitle: "SynBio Ready-to-Order CDS FASTA (100% Identical Protein)",
    copyFastaBtn: "Copy FASTA",
    loadHarmBtn: "Load Designed Gene into Simulator",
    phaseTitle: "TASEP Statistical-Physics Phase Transition Diagram (α Sweep)",
    phaseSub: "Sweeps initiation rate α from 0.08 to 2.20 /s to reveal the exact critical threshold α_c where the transcript transitions from Initiation-Limited (Low-Density) into a Ribosome Traffic Jam (High-Density).",
    runPhaseBtn: "📈 Run Phase Transition Sweep",
    legFluxWT: "Protein Flux J(α) — Current Gene",
    legNgdWT: "NGD Collision Risk — Current Gene",
    legFluxOpt: "Flux J(α) — SynBio Optimized (if run)",
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
  lastMeltResult: null,
  lastSynBioPayload: null,
  lastPhaseData: null,
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

  const sampleSel = $("sample-select");
  state.meta.samples.forEach((s) => {
    const opt = document.createElement("option");
    opt.value = s.id;
    opt.textContent = s.name;
    sampleSel.appendChild(opt);
  });

  ["organism-select", "source-org-select", "target-org-select"].forEach((id) => {
    const sel = $(id);
    state.meta.organisms.forEach((o) => {
      const opt = document.createElement("option");
      opt.value = o.key;
      opt.textContent = `${o.name} (${o.common_name})`;
      sel.appendChild(opt);
    });
  });

  sampleSel.value = "polybasic_arrest";
  $("organism-select").value = "ecoli";
  $("source-org-select").value = "human";
  $("target-org-select").value = "ecoli";

  const first = state.meta.samples.find((s) => s.id === "polybasic_arrest") || state.meta.samples[0];
  state.currentSequence = first.sequence;
  state.currentName = first.name;

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
    state.lastMeltResult = null;
    state.lastSynBioPayload = null;
    $("organism-select").value = s.recommended_host;
    $("source-org-select").value = s.organism;
    $("target-org-select").value = s.recommended_host;
    $("opt-vaccine").checked = s.id === "sars_cov2_rbd";
    const hostObj = state.meta.organisms.find((o) => o.key === s.recommended_host);
    if (hostObj) {
      state.alpha = hostObj.default_initiation_rate;
      $("alpha-slider").value = state.alpha;
      $("alpha-value").textContent = `${state.alpha.toFixed(2)} /s`;
    }
    updateSampleDescription();
    $("mutation-result-card").classList.add("hidden");
    $("hotspots-panel").classList.add("hidden");
    $("melt-comparison").classList.add("hidden");
    $("synbio-results").classList.add("hidden");
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
    await runSimulation();
  });

  $("alpha-slider").addEventListener("input", () => {
    state.alpha = parseFloat($("alpha-slider").value);
    $("alpha-value").textContent = `${state.alpha.toFixed(2)} /s`;
  });
  $("alpha-slider").addEventListener("change", () => runSimulation());
  $("toggle-5prime").addEventListener("change", () => runSimulation());
  $("toggle-tunnel").addEventListener("change", () => runSimulation());
  $("run-btn").addEventListener("click", () => runSimulation());

  // Live NCBI RefSeq Fetcher (CORS-enabled NCBI E-utilities!)
  $("toggle-ncbi-btn").addEventListener("click", () => {
    $("ncbi-box").classList.toggle("hidden");
    $("custom-fasta-box").classList.add("hidden");
  });
  $("ncbi-fetch-btn").addEventListener("click", fetchFromNCBI);

  $("toggle-fasta-btn").addEventListener("click", () => {
    $("custom-fasta-box").classList.toggle("hidden");
    $("ncbi-box").classList.add("hidden");
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
        drawTunnelChart();
      } else if (target === "folding" && state.simulation) {
        renderFoldingTab(state.simulation.folding_5prime);
      }
    });
  });

  // Playback
  $("play-pause-btn").addEventListener("click", () => {
    state.playing = !state.playing;
    $("play-pause-btn").textContent = state.playing ? I18N[state.lang].pause : I18N[state.lang].play;
  });
  $("speed-select").addEventListener("change", () => { state.speed = parseFloat($("speed-select").value); });

  // PNG Figure Exporters
  document.querySelectorAll(".png-btn").forEach((btn) => {
    btn.addEventListener("click", () => {
      const c = $(btn.dataset.canvas);
      if (!c) return;
      const a = document.createElement("a");
      a.download = `${btn.dataset.name || "riboflow_figure"}.png`;
      a.href = c.toDataURL("image/png");
      a.click();
    });
  });

  // Printable Lab Report
  $("export-report-btn").addEventListener("click", exportPrintableReport);

  // 5' Hairpin Melter
  $("melt-hairpin-btn").addEventListener("click", runHairpinMelt);
  $("load-melted-btn").addEventListener("click", async () => {
    if (!state.lastMeltResult) return;
    state.currentSequence = state.lastMeltResult.unfolded_sequence;
    state.currentName = `${state.currentName} (5'-Melted)`;
    $("melt-comparison").classList.add("hidden");
    document.querySelector('[data-tab="traffic"]').click();
    await runSimulation();
  });

  // Silent Mutation Lab
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

  // SynBio & mRNA Vaccine Designer
  $("run-synbio-btn").addEventListener("click", runSynBioDesign);
  $("copy-harm-btn").addEventListener("click", async () => {
    if (!state.lastSynBioPayload) return;
    await navigator.clipboard.writeText($("harm-fasta-out").textContent);
    const btn = $("copy-harm-btn");
    const prev = btn.textContent;
    btn.textContent = I18N[state.lang].copied;
    setTimeout(() => { btn.textContent = prev; }, 1400);
  });
  $("load-harm-btn").addEventListener("click", async () => {
    if (!state.lastSynBioPayload) return;
    const sb = state.lastSynBioPayload.synbio;
    state.currentSequence = sb.designed_sequence;
    state.currentName = `${state.currentName} (SynBio v2)`;
    $("organism-select").value = sb.host_organism;
    state.organism = sb.host_organism;
    state.edits = {};
    document.querySelector('[data-tab="traffic"]').click();
    await runSimulation();
  });

  // Phase Diagram Sweep
  $("run-phase-btn").addEventListener("click", runPhaseSweep);
}

/* ---------- LIVE NCBI REFSEQ FETCHER ---------- */
function extractBestCDS(rawFasta) {
  const clean = rawFasta
    .split(/\r?\n/)
    .filter((l) => !l.startsWith(">"))
    .join("")
    .toUpperCase()
    .replace(/[^ATGC]/g, "");
  // Scan all 3 frames for the longest valid ATG...STOP ORF (capped to 450 codons for fast browser simulation)
  let best = "";
  for (let frame = 0; frame < 3; frame++) {
    let i = frame;
    while (i + 3 <= clean.length) {
      if (clean.slice(i, i + 3) === "ATG") {
        for (let j = i + 3; j + 3 <= clean.length; j += 3) {
          const c = clean.slice(j, j + 3);
          if (c === "TAA" || c === "TAG" || c === "TGA") {
            const orf = clean.slice(i, j + 3);
            if (orf.length > best.length && orf.length >= 45) best = orf;
            i = j + 3;
            break;
          }
        }
      }
      i += 3;
    }
  }
  if (!best) throw new Error("No valid coding ORF (ATG...STOP) found in NCBI record.");
  // Cap very long transcripts at 350 codons + stop for interactive responsiveness
  if (best.length > 1053) {
    best = best.slice(0, 1050) + "TAA";
  }
  return best;
}

async function fetchFromNCBI() {
  const q = $("ncbi-query").value.trim();
  const statusEl = $("ncbi-status");
  if (!q) return;
  $("ncbi-fetch-btn").disabled = true;
  statusEl.textContent = "Querying NCBI E-utilities...";
  try {
    const term = /^(NM_|XM_|NC_)/i.test(q)
      ? q
      : `${q}[Gene] AND Homo sapiens[Organism] AND RefSeq[Filter] AND mRNA[Filter]`;
    const sRes = await fetch(
      `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=${encodeURIComponent(term)}&retmax=1&retmode=json`
    );
    const sData = await sRes.json();
    const id = sData.esearchresult && sData.esearchresult.idlist && sData.esearchresult.idlist[0];
    if (!id) throw new Error("Gene not found in NCBI RefSeq.");

    statusEl.textContent = `Downloading FASTA (ID ${id})...`;
    const fRes = await fetch(
      `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${id}&rettype=fasta&retmode=text`
    );
    const fastaText = await fRes.text();
    const cds = extractBestCDS(fastaText);

    state.currentSequence = cds;
    state.currentName = `NCBI:${q.toUpperCase()} (${cds.length / 3} codons)`;
    state.edits = {};
    statusEl.textContent = `✓ Loaded ${cds.length / 3} codons from NCBI!`;
    await runSimulation();
  } catch (err) {
    statusEl.textContent = `Error: ${err.message}`;
  } finally {
    $("ncbi-fetch-btn").disabled = false;
  }
}

/* ---------- CORE SIMULATION ---------- */
async function runSimulation() {
  state.organism = $("organism-select").value;
  state.alpha = parseFloat($("alpha-slider").value);

  const data = await apiPost("/api/simulate", {
    sequence: state.currentSequence,
    organism: state.organism,
    initiation_rate: state.alpha,
    duration: 110,
    couple_5prime_folding: $("toggle-5prime").checked,
    enable_tunnel_drag: $("toggle-tunnel").checked,
    seed: 42,
  });

  state.analysis = data.analysis;
  state.simulation = data.simulation;
  state.currentSequence = data.analysis.sequence;
  state.frameCursor = 0;

  $("active-gene-pill").textContent = `${state.currentName} • ${data.analysis.codon_count} codons`;

  updateKPIs();
  drawRiboseqChart();
  drawMinmaxChart();
  drawTunnelChart();
  renderFoldingTab(state.simulation.folding_5prime);
  renderCodonStrip();
}

function updateKPIs() {
  const sim = state.simulation;
  const an = state.analysis;
  const f5 = sim.folding_5prime;

  $("kpi-flux").textContent = sim.protein_flux_per_min.toFixed(1);
  $("kpi-eff-alpha").textContent = `α_eff = ${sim.effective_initiation_rate.toFixed(2)} /s`;

  const mfeEl = $("kpi-mfe");
  mfeEl.textContent = `${f5.mfe_kcal_mol.toFixed(1)} kcal`;
  mfeEl.className = "kpi-val " + (f5.initiation_accessibility < 0.5 ? "risk-high" : f5.initiation_accessibility < 0.8 ? "risk-mid" : "risk-low");
  $("kpi-access").textContent = `Access ${(f5.initiation_accessibility * 100).toFixed(0)}% • ${f5.arcs.length} bp`;

  $("kpi-collision").textContent = `${(sim.collision_fraction * 100).toFixed(1)}%`;
  $("kpi-collision-events").textContent = `${sim.collision_events} events • peak #${sim.bottleneck_codon + 1}`;

  const ngdEl = $("kpi-ngd");
  ngdEl.textContent = `${sim.ngd_risk_score.toFixed(1)} / 100`;
  ngdEl.className = "kpi-val " + (sim.ngd_risk_score >= 45 ? "risk-high" : sim.ngd_risk_score >= 18 ? "risk-mid" : "risk-low");
  $("kpi-queue").textContent = `max queue: ${sim.max_queue_length} ribos`;

  const tunEl = $("kpi-tunnel");
  tunEl.textContent = `+${sim.peak_tunnel_charge.toFixed(1)} q`;
  tunEl.className = "kpi-val " + (sim.peak_tunnel_charge >= 6 ? "risk-high" : sim.peak_tunnel_charge >= 3.5 ? "risk-mid" : "risk-low");
  $("kpi-cai-gc").textContent = `CAI ${an.cai.toFixed(3)} • GC3 ${an.gc3_percent}%`;

  $("kpi-phase").textContent = sim.phase_regime;
  $("kpi-ribos").textContent = `avg ${sim.mean_ribosomes_on_mrna} ribos on mRNA`;
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
      if (state.frameCursor >= state.simulation.frames.length - 1) state.frameCursor = 0;
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

  ctx.strokeStyle = "rgba(35, 50, 77, 0.35)";
  ctx.lineWidth = 1;
  for (let x = 50; x < W - 40; x += 80) {
    ctx.beginPath(); ctx.moveTo(x, 20); ctx.lineTo(x, H - 25); ctx.stroke();
  }

  const padL = 54, padR = 54;
  const trackW = W - padL - padR;
  const trackY = 178;
  const L = sim.codon_count;
  const codonW = trackW / L;
  const maxRate = Math.max(...sim.codon_rates);

  // Shade polybasic exit-tunnel drag zones along track
  if (sim.enable_tunnel_drag && sim.tunnel_drag_factors) {
    for (let i = 0; i < L; i++) {
      if (sim.tunnel_drag_factors[i] < 0.78) {
        const x = padL + i * codonW;
        ctx.fillStyle = "rgba(167, 139, 250, 0.22)";
        ctx.fillRect(x, trackY - 34, Math.max(2, codonW), 68);
      }
    }
  }

  ctx.fillStyle = "#8ea3c4";
  ctx.font = "bold 12px monospace";
  ctx.fillText("5' AUG", 8, trackY + 4);
  ctx.fillText("3' STOP", W - 50, trackY + 4);

  for (let i = 0; i < L; i++) {
    const x = padL + i * codonW;
    ctx.fillStyle = codonSpeedColor(sim.codon_rates[i], maxRate);
    ctx.fillRect(x, trackY - 5, Math.max(1.5, codonW - 0.6), 10);
  }

  const bx = padL + (sim.bottleneck_codon + 0.5) * codonW;
  ctx.strokeStyle = "rgba(244, 63, 94, 0.55)";
  ctx.setLineDash([4, 4]);
  ctx.beginPath(); ctx.moveTo(bx, 28); ctx.lineTo(bx, trackY + 38); ctx.stroke();
  ctx.setLineDash([]);
  ctx.fillStyle = "#fda4af";
  ctx.font = "11px monospace";
  ctx.fillText(`Bottleneck #${sim.bottleneck_codon + 1}`, Math.min(W - 140, Math.max(padL, bx - 42)), trackY + 54);

  const nextMap = new Map(frameB.ribosomes.map((r) => [r.id, r]));

  frameA.ribosomes.forEach((r) => {
    const next = nextMap.get(r.id);
    const interpPos = next ? r.pos + (next.pos - r.pos) * frac : r.pos;
    const cx = padL + (interpPos + 0.5) * codonW;
    const colliding = r.colliding;
    const tunnelStall = r.tunnel_stalled;

    const fpW = Math.max(18, sim.ribosome_footprint * codonW);
    ctx.fillStyle = colliding
      ? "rgba(244, 63, 94, 0.22)"
      : tunnelStall
      ? "rgba(167, 139, 250, 0.22)"
      : "rgba(45, 212, 191, 0.12)";
    ctx.fillRect(cx - fpW * 0.5, trackY - 26, fpW, 52);

    // Growing nascent polypeptide chain
    const progress = interpPos / Math.max(1, L - 1);
    const beads = Math.max(2, Math.floor(progress * 11));
    ctx.strokeStyle = colliding ? "#fda4af" : tunnelStall ? "#c4b5fd" : "#7dd3fc";
    ctx.lineWidth = tunnelStall ? 2.8 : 2;
    ctx.beginPath();
    let px = cx, py = trackY - 26;
    ctx.moveTo(px, py);
    for (let b = 1; b <= beads; b++) {
      px += Math.sin(b * 1.3 + r.id) * 6;
      py -= 9;
      ctx.lineTo(px, py);
    }
    ctx.stroke();

    // Ribosome subunits
    ctx.fillStyle = colliding ? "#f43f5e" : tunnelStall ? "#a78bfa" : "#2dd4bf";
    ctx.beginPath();
    ctx.ellipse(cx, trackY - 15, Math.min(20, fpW * 0.44), 12, 0, 0, Math.PI * 2);
    ctx.fill();

    ctx.fillStyle = colliding ? "#9f1239" : tunnelStall ? "#5b21b6" : "#0f766e";
    ctx.beginPath();
    ctx.ellipse(cx, trackY + 11, Math.min(15, fpW * 0.34), 7, 0, 0, Math.PI * 2);
    ctx.fill();

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
  const W = canvas.width, H = canvas.height;
  ctx.clearRect(0, 0, W, H);
  const pad = { l: 38, r: 12, t: 16, b: 24 };
  const plotW = W - pad.l - pad.r, plotH = H - pad.t - pad.b;
  const L = sim.codon_occupancy.length;
  const maxY = Math.max(0.12, ...sim.codon_occupancy, ...sim.collision_profile);

  ctx.strokeStyle = "#23324d";
  ctx.beginPath(); ctx.moveTo(pad.l, pad.t); ctx.lineTo(pad.l, H - pad.b); ctx.lineTo(W - pad.r, H - pad.b); ctx.stroke();

  const fillSeries = (arr, fill, stroke) => {
    ctx.fillStyle = fill; ctx.strokeStyle = stroke; ctx.lineWidth = 1.7;
    ctx.beginPath(); ctx.moveTo(pad.l, H - pad.b);
    for (let i = 0; i < L; i++) {
      const x = pad.l + (i / Math.max(1, L - 1)) * plotW;
      const y = H - pad.b - (arr[i] / maxY) * plotH;
      ctx.lineTo(x, y);
    }
    ctx.lineTo(W - pad.r, H - pad.b); ctx.closePath(); ctx.fill(); ctx.stroke();
  };
  fillSeries(sim.codon_occupancy, "rgba(45, 212, 191, 0.25)", "#2dd4bf");
  fillSeries(sim.collision_profile, "rgba(244, 63, 94, 0.48)", "#f43f5e");

  ctx.fillStyle = "#8ea3c4"; ctx.font = "10px monospace";
  ctx.fillText(maxY.toFixed(2), 4, pad.t + 8);
  ctx.fillText("1", pad.l, H - 6); ctx.fillText(`${L}`, W - 26, H - 6);
}

function drawMinmaxChart() {
  const canvas = $("minmax-canvas");
  const ctx = canvas.getContext("2d");
  const an = state.analysis;
  if (!an) return;
  const W = canvas.width, H = canvas.height;
  ctx.clearRect(0, 0, W, H);
  const pad = { l: 38, r: 12, t: 16, b: 24 };
  const plotW = W - pad.l - pad.r, plotH = H - pad.t - pad.b;
  const midY = pad.t + plotH / 2;
  const scores = an.codons.map((c) => c.minmax_score);
  const L = scores.length;

  ctx.strokeStyle = "#334155"; ctx.setLineDash([4, 4]);
  ctx.beginPath(); ctx.moveTo(pad.l, midY); ctx.lineTo(W - pad.r, midY); ctx.stroke(); ctx.setLineDash([]);

  for (let i = 0; i < L - 1; i++) {
    const x1 = pad.l + (i / (L - 1)) * plotW, y1 = midY - (scores[i] / 100) * (plotH / 2);
    const x2 = pad.l + ((i + 1) / (L - 1)) * plotW, y2 = midY - (scores[i + 1] / 100) * (plotH / 2);
    ctx.strokeStyle = scores[i] < -10 ? "#f59e0b" : "#38bdf8";
    ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x2, y2); ctx.stroke();
  }
  ctx.fillStyle = "#8ea3c4"; ctx.font = "10px monospace";
  ctx.fillText("+100", 4, pad.t + 8); ctx.fillText("   0", 4, midY + 4); ctx.fillText("-100", 4, H - pad.b);
}

function drawTunnelChart() {
  const canvas = $("tunnel-canvas");
  const ctx = canvas.getContext("2d");
  const sim = state.simulation;
  if (!sim || !sim.tunnel_charge_profile) return;
  const W = canvas.width, H = canvas.height;
  ctx.clearRect(0, 0, W, H);
  const pad = { l: 36, r: 12, t: 16, b: 24 };
  const plotW = W - pad.l - pad.r, plotH = H - pad.t - pad.b;
  const arr = sim.tunnel_charge_profile;
  const L = arr.length;
  const maxY = Math.max(6, ...arr);

  ctx.strokeStyle = "#23324d";
  ctx.beginPath(); ctx.moveTo(pad.l, pad.t); ctx.lineTo(pad.l, H - pad.b); ctx.lineTo(W - pad.r, H - pad.b); ctx.stroke();

  // Drag threshold line at +3.0 charge
  const threshY = H - pad.b - (3.0 / maxY) * plotH;
  ctx.strokeStyle = "rgba(244, 63, 94, 0.55)"; ctx.setLineDash([4, 4]);
  ctx.beginPath(); ctx.moveTo(pad.l, threshY); ctx.lineTo(W - pad.r, threshY); ctx.stroke(); ctx.setLineDash([]);

  ctx.fillStyle = "rgba(167, 139, 250, 0.32)"; ctx.strokeStyle = "#a78bfa"; ctx.lineWidth = 1.9;
  ctx.beginPath(); ctx.moveTo(pad.l, H - pad.b);
  for (let i = 0; i < L; i++) {
    const x = pad.l + (i / Math.max(1, L - 1)) * plotW;
    const y = H - pad.b - (arr[i] / maxY) * plotH;
    ctx.lineTo(x, y);
  }
  ctx.lineTo(W - pad.r, H - pad.b); ctx.closePath(); ctx.fill(); ctx.stroke();

  ctx.fillStyle = "#8ea3c4"; ctx.font = "10px monospace";
  ctx.fillText(`+${maxY.toFixed(0)}`, 4, pad.t + 8);
  ctx.fillText("+3q", 4, threshY + 4);
  ctx.fillText("Drag threshold (+3 K/R)", pad.l + 6, threshY - 5);
}

/* ---------- TAB 2: 5' mRNA HAIRPIN THERMODYNAMICS ---------- */
function renderFoldingTab(f5) {
  if (!f5) return;
  $("fold-metrics").innerHTML = `
    <div class="mut-stat"><span>5' Window MFE (ΔG)</span><strong class="${f5.initiation_accessibility < 0.5 ? "risk-high" : "risk-low"}">${f5.mfe_kcal_mol.toFixed(2)} kcal/mol</strong></div>
    <div class="mut-stat"><span>Initiation Accessibility</span><strong>${(f5.initiation_accessibility * 100).toFixed(1)}%</strong></div>
    <div class="mut-stat"><span>AUG Zone (nt 1-15) Paired</span><strong>${(f5.start_region_paired * 100).toFixed(1)}%</strong></div>
    <div class="mut-stat"><span>Verdict</span><strong>${f5.severity}</strong></div>
  `;
  drawArcDiagram($("arc-canvas"), f5);
  $("dot-bracket-box").textContent = `5' RNA (${f5.window_nt} nt): ${f5.rna_sequence}\nDot-Bracket  : ${f5.dot_bracket}`;
}

function drawArcDiagram(canvas, f5) {
  const ctx = canvas.getContext("2d");
  const W = canvas.width, H = canvas.height;
  ctx.clearRect(0, 0, W, H);
  const padL = 32, padR = 32;
  const n = f5.rna_sequence.length;
  const step = (W - padL - padR) / Math.max(1, n - 1);
  const baseY = H - 42;

  // Shade AUG initiation zone (first 15 nt)
  ctx.fillStyle = "rgba(244, 63, 94, 0.12)";
  ctx.fillRect(padL - 8, 12, step * Math.min(15, n) + 8, baseY - 4);
  ctx.fillStyle = "#fda4af"; ctx.font = "11px monospace";
  ctx.fillText("AUG Zone (nt 1-15)", padL, 26);

  // Draw base-pair arcs
  f5.arcs.forEach((arc) => {
    const x1 = padL + arc.i * step;
    const x2 = padL + arc.j * step;
    const cx = (x1 + x2) / 2;
    const rx = (x2 - x1) / 2;
    const ry = Math.min(baseY - 34, rx * 0.72 + 12);
    ctx.strokeStyle =
      arc.pair === "G-C" || arc.pair === "C-G"
        ? "#10b981"
        : arc.pair === "A-U" || arc.pair === "U-A"
        ? "#38bdf8"
        : "#f59e0b";
    ctx.lineWidth = 2.2;
    ctx.beginPath();
    ctx.ellipse(cx, baseY - 12, rx, ry, 0, Math.PI, 0, false);
    ctx.stroke();
  });

  // Draw RNA sequence bases
  const baseColor = { A: "#38bdf8", U: "#f43f5e", G: "#10b981", C: "#f59e0b" };
  ctx.font = "bold 13px monospace";
  ctx.textAlign = "center";
  for (let i = 0; i < n; i++) {
    const ch = f5.rna_sequence[i];
    const x = padL + i * step;
    ctx.fillStyle = baseColor[ch] || "#e2e8f0";
    ctx.fillText(ch, x, baseY + 6);
    if (i % 5 === 0) {
      ctx.fillStyle = "#64748b";
      ctx.font = "10px monospace";
      ctx.fillText(`${i + 1}`, x, baseY + 24);
      ctx.font = "bold 13px monospace";
    }
  }
  ctx.textAlign = "left";
}

async function runHairpinMelt() {
  const btn = $("melt-hairpin-btn");
  btn.disabled = true;
  try {
    const res = await apiPost("/api/unfold", { sequence: state.currentSequence });
    state.lastMeltResult = res;
    $("melt-comparison").classList.remove("hidden");
    $("melt-gain-badge").textContent =
      `✓ ΔMFE ${res.delta_mfe >= 0 ? "+" : ""}${res.delta_mfe} kcal/mol • Access ${(res.before.initiation_accessibility * 100).toFixed(0)}% → ${(res.after.initiation_accessibility * 100).toFixed(0)}% (${res.accessibility_gain_percent >= 0 ? "+" : ""}${res.accessibility_gain_percent}%) • ${res.codons_changed} codons changed`;
    drawArcDiagram($("arc-after-canvas"), res.after);
    $("dot-bracket-after").textContent =
      `5' RNA (${res.after.window_nt} nt): ${res.after.rna_sequence}\nDot-Bracket  : ${res.after.dot_bracket}`;
  } finally {
    btn.disabled = false;
  }
}

/* ---------- TAB 3: SILENT MUTATION LAB ---------- */
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
    if (!locked) btn.addEventListener("click", () => openCodonEditor(c));
    strip.appendChild(btn);
  });
}

function openCodonEditor(codonObj) {
  state.selectedCodon = codonObj.index;
  renderCodonStrip();
  $("codon-editor-drawer").classList.remove("hidden");
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
      if (syn === codonObj.codon) delete state.edits[codonObj.index];
      else state.edits[codonObj.index] = syn;
      renderCodonStrip();
      openCodonEditor(codonObj);
      if (Object.keys(state.edits).length > 0) await runMutationComparison();
      else $("mutation-result-card").classList.add("hidden");
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
  const fluxSign = comp.delta_flux_per_min >= 0 ? "+" : "";
  const ngdSign = comp.delta_ngd_risk >= 0 ? "+" : "";
  const collDeltaPct = (comp.delta_collision_fraction * 100).toFixed(2);

  $("mut-metrics").innerHTML = `
    <div class="mut-stat"><span>Protein Flux (J)</span><strong>${comp.wildtype_sim.protein_flux_per_min} → ${comp.mutant_sim.protein_flux_per_min} (${fluxSign}${comp.delta_flux_percent}%)</strong></div>
    <div class="mut-stat"><span>Disome Collision %</span><strong>${(comp.wildtype_sim.collision_fraction * 100).toFixed(1)}% → ${(comp.mutant_sim.collision_fraction * 100).toFixed(1)}% (${collDeltaPct >= 0 ? "+" : ""}${collDeltaPct}%)</strong></div>
    <div class="mut-stat"><span>NGD / RQC Risk</span><strong class="${comp.delta_ngd_risk > 5 ? "risk-high" : "risk-low"}">${comp.wildtype_sim.ngd_risk_score} → ${comp.mutant_sim.ngd_risk_score} (${ngdSign}${comp.delta_ngd_risk})</strong></div>
    <div class="mut-stat"><span>CAI Change</span><strong>${comp.wildtype_cai.toFixed(3)} → ${comp.mutant_cai.toFixed(3)}</strong></div>
  `;
  drawMutationDiffCanvas(comp);
}

function drawMutationDiffCanvas(comp) {
  const canvas = $("mut-diff-canvas");
  const ctx = canvas.getContext("2d");
  const W = canvas.width, H = canvas.height;
  ctx.clearRect(0, 0, W, H);
  const pad = { l: 48, r: 20, t: 22, b: 26 };
  const plotW = W - pad.l - pad.r, plotH = H - pad.t - pad.b;
  const wt = comp.wildtype_sim.codon_occupancy, mut = comp.mutant_sim.codon_occupancy;
  const L = wt.length, maxY = Math.max(0.12, ...wt, ...mut);

  ctx.strokeStyle = "#23324d";
  ctx.beginPath(); ctx.moveTo(pad.l, pad.t); ctx.lineTo(pad.l, H - pad.b); ctx.lineTo(W - pad.r, H - pad.b); ctx.stroke();

  const drawLine = (arr, color, dashed) => {
    ctx.strokeStyle = color; ctx.lineWidth = 2; ctx.setLineDash(dashed ? [5, 4] : []);
    ctx.beginPath();
    arr.forEach((v, i) => {
      const x = pad.l + (i / (L - 1)) * plotW, y = H - pad.b - (v / maxY) * plotH;
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    });
    ctx.stroke(); ctx.setLineDash([]);
  };
  drawLine(wt, "#94a3b8", true); drawLine(mut, "#f43f5e", false);
  comp.edits.forEach((e) => {
    const x = pad.l + (e.index / (L - 1)) * plotW;
    ctx.strokeStyle = "#f59e0b"; ctx.beginPath(); ctx.moveTo(x, pad.t); ctx.lineTo(x, H - pad.b); ctx.stroke();
    ctx.fillStyle = "#fbbf24"; ctx.font = "11px monospace";
    ctx.fillText(`#${e.index + 1} ${e.from_codon}→${e.to_codon}`, Math.min(W - 140, x + 4), pad.t + 12);
  });
}

async function scanHotspots() {
  const btn = $("scan-hotspots-btn");
  btn.disabled = true;
  try {
    const data = await apiPost("/api/scan-mutations", {
      sequence: state.currentSequence, organism: state.organism, initiation_rate: state.alpha, top_k: 6, seed: 42,
    });
    const list = $("hotspots-list");
    list.innerHTML = "";
    $("hotspots-panel").classList.remove("hidden");
    data.hotspots.forEach((h) => {
      const el = document.createElement("div");
      el.className = "hotspot-item";
      el.innerHTML = `<strong>#${h.index + 1} ${h.aa_three}: ${h.from_codon} → ${h.to_codon}</strong><div class="muted">Rate -${h.rate_drop_percent}% • ΔNGD <span class="risk-high">${h.delta_ngd_risk >= 0 ? "+" : ""}${h.delta_ngd_risk}</span> • ΔFlux ${h.delta_flux_per_min >= 0 ? "+" : ""}${h.delta_flux_per_min}</div>`;
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

/* ---------- TAB 4: MULTI-OBJECTIVE SYNBIO & mRNA VACCINE DESIGNER ---------- */
async function runSynBioDesign() {
  const btn = $("run-synbio-btn");
  btn.disabled = true;
  try {
    const payload = await apiPost("/api/synbio", {
      sequence: state.currentSequence,
      source: $("source-org-select").value,
      host: $("target-org-select").value,
      preserve_pauses: $("opt-pauses").checked,
      unfold_5prime: $("opt-unfold").checked,
      uridine_deplete: $("opt-vaccine").checked,
      scrub_restriction_sites: $("opt-scrub").checked,
      avoid_homopolymers: $("opt-hp").checked,
      initiation_rate: state.alpha,
      seed: 42,
    });
    state.lastSynBioPayload = payload;
    const sb = payload.synbio;
    const hr = payload.harmonization;

    $("synbio-results").classList.remove("hidden");
    $("synbio-scorecard").innerHTML = `
      <div class="mut-stat">
        <span>5' Hairpin MFE (ΔG)</span>
        <strong class="risk-low">${sb.before_qc.folding_5prime.mfe_kcal_mol.toFixed(1)} → ${sb.after_qc.folding_5prime.mfe_kcal_mol.toFixed(1)} kcal/mol</strong>
      </div>
      <div class="mut-stat">
        <span>Uridine / UU (TLR7-8)</span>
        <strong class="risk-low">${sb.before_qc.uridine_percent}% (${sb.before_qc.uu_dinucleotides} UU) → ${sb.after_qc.uridine_percent}% (${sb.after_qc.uu_dinucleotides} UU)</strong>
      </div>
      <div class="mut-stat">
        <span>Cloning Restriction Sites</span>
        <strong class="${sb.after_qc.restriction_hits.length === 0 ? "risk-low" : "risk-mid"}">${sb.before_qc.restriction_hits.length} → ${sb.after_qc.restriction_hits.length} sites</strong>
      </div>
      <div class="mut-stat">
        <span>Homopolymers (≥5 nt)</span>
        <strong class="${sb.after_qc.homopolymer_violations === 0 ? "risk-low" : "risk-mid"}">${sb.before_qc.homopolymer_violations} → ${sb.after_qc.homopolymer_violations}</strong>
      </div>
    `;

    const tbody = $("synbio-tbody");
    tbody.innerHTML = "";
    const rows = [
      {
        label: hr.wildtype_in_host.label,
        cai: hr.wildtype_in_host.cai_in_host,
        r: hr.wildtype_in_host.rhythm_correlation,
        mfe: sb.before_qc.folding_5prime.mfe_kcal_mol,
        sim: sb.before_sim,
        hi: false,
      },
      {
        label: hr.naive_optimized.label,
        cai: hr.naive_optimized.cai_in_host,
        r: hr.naive_optimized.rhythm_correlation,
        mfe: hr.naive_optimized.simulation.folding_5prime.mfe_kcal_mol,
        sim: hr.naive_optimized.simulation,
        hi: false,
      },
      {
        label: "RiboFlow v2.0 SynBio Designed ★",
        cai: sb.designed_cai,
        r: sb.designed_rhythm_r,
        mfe: sb.after_qc.folding_5prime.mfe_kcal_mol,
        sim: sb.after_sim,
        hi: true,
      },
    ];

    rows.forEach((row) => {
      const tr = document.createElement("tr");
      if (row.hi) tr.className = "highlight-row";
      tr.innerHTML = `
        <td><strong>${row.label}</strong></td>
        <td>${row.cai.toFixed(3)}</td>
        <td>${row.r >= 0 ? "+" : ""}${row.r.toFixed(3)}</td>
        <td>${row.mfe.toFixed(1)}</td>
        <td>${row.sim.protein_flux_per_min.toFixed(2)}</td>
        <td>${(row.sim.collision_fraction * 100).toFixed(2)}%</td>
        <td>${row.sim.ngd_risk_score.toFixed(1)}</td>
      `;
      tbody.appendChild(tr);
    });

    drawHarmonizationCanvas(sb.native_source_minmax, hr.naive_optimized.minmax_profile, sb.after_minmax);
    $("harm-fasta-out").textContent =
      `>${state.currentName}_riboflow_v2_${sb.source_organism}_to_${sb.host_organism}\n${sb.designed_sequence}`;
  } finally {
    btn.disabled = false;
  }
}

function drawHarmonizationCanvas(nativeArr, naiveArr, synArr) {
  const canvas = $("harm-canvas");
  const ctx = canvas.getContext("2d");
  const W = canvas.width, H = canvas.height;
  ctx.clearRect(0, 0, W, H);
  const pad = { l: 46, r: 18, t: 18, b: 26 };
  const plotW = W - pad.l - pad.r, plotH = H - pad.t - pad.b;
  const midY = pad.t + plotH / 2;

  ctx.strokeStyle = "#334155"; ctx.setLineDash([4, 4]);
  ctx.beginPath(); ctx.moveTo(pad.l, midY); ctx.lineTo(W - pad.r, midY); ctx.stroke(); ctx.setLineDash([]);

  const drawSeries = (arr, color, width, dashed) => {
    ctx.strokeStyle = color; ctx.lineWidth = width; ctx.setLineDash(dashed ? [6, 4] : []);
    ctx.beginPath();
    arr.forEach((v, i) => {
      const x = pad.l + (i / Math.max(1, arr.length - 1)) * plotW;
      const y = midY - (v / 100) * (plotH / 2);
      if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    });
    ctx.stroke(); ctx.setLineDash([]);
  };

  drawSeries(naiveArr, "#a78bfa", 1.8, false);
  drawSeries(nativeArr, "#e2e8f0", 2.2, true);
  drawSeries(synArr, "#2dd4bf", 2.5, false);

  ctx.fillStyle = "#8ea3c4"; ctx.font = "11px monospace";
  ctx.fillText("+100", 6, pad.t + 8); ctx.fillText("   0", 6, midY + 4); ctx.fillText("-100", 6, H - pad.b);
}

/* ---------- TAB 5: TASEP PHASE TRANSITION DIAGRAM ---------- */
async function runPhaseSweep() {
  const btn = $("run-phase-btn");
  btn.disabled = true;
  try {
    const req = {
      sequence: state.currentSequence,
      organism: state.organism,
      couple_5prime_folding: $("toggle-5prime").checked,
      enable_tunnel_drag: $("toggle-tunnel").checked,
      seed: 42,
    };
    if (state.lastSynBioPayload && state.lastSynBioPayload.synbio.designed_sequence !== state.currentSequence) {
      req.compare_sequence = state.lastSynBioPayload.synbio.designed_sequence;
    }
    const data = await apiPost("/api/phase", req);
    state.lastPhaseData = data;
    const wt = data.wildtype;
    const opt = data.optimized;

    $("phase-metrics").classList.remove("hidden");
    $("phase-metrics").innerHTML = `
      <div class="mut-stat"><span>Critical Jamming α_c (Current)</span><strong class="${wt.critical_alpha && wt.critical_alpha < 0.5 ? "risk-high" : "risk-low"}">${wt.critical_alpha ? wt.critical_alpha + " /s" : "Jam-Free"}</strong></div>
      <div class="mut-stat"><span>Peak Protein Flux (Current)</span><strong>${wt.max_flux_per_min.toFixed(1)} /min</strong></div>
      <div class="mut-stat"><span>SynBio Comparison</span><strong>${opt ? `Peak ${opt.max_flux_per_min.toFixed(1)} /min` : "Run Tab 4 to overlay"}</strong></div>
      <div class="mut-stat"><span>Host Organism</span><strong>${wt.organism.toUpperCase()}</strong></div>
    `;
    drawPhaseCanvas(wt, opt);
  } finally {
    btn.disabled = false;
  }
}

function drawPhaseCanvas(wt, opt) {
  const canvas = $("phase-canvas");
  const ctx = canvas.getContext("2d");
  const W = canvas.width, H = canvas.height;
  ctx.clearRect(0, 0, W, H);
  const pad = { l: 58, r: 58, t: 24, b: 36 };
  const plotW = W - pad.l - pad.r, plotH = H - pad.t - pad.b;

  const maxAlpha = Math.max(...wt.points.map((p) => p.alpha));
  const maxFlux = Math.max(
    10,
    ...wt.points.map((p) => p.protein_flux_per_min),
    ...(opt ? opt.points.map((p) => p.protein_flux_per_min) : [])
  ) * 1.15;

  const xOf = (a) => pad.l + (a / maxAlpha) * plotW;
  const yFlux = (j) => H - pad.b - (j / maxFlux) * plotH;
  const yNgd = (n) => H - pad.b - (n / 100) * plotH;

  // Grid & axes
  ctx.strokeStyle = "#23324d"; ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(pad.l, pad.t); ctx.lineTo(pad.l, H - pad.b); ctx.lineTo(W - pad.r, H - pad.b); ctx.lineTo(W - pad.r, pad.t);
  ctx.stroke();

  // Critical alpha_c vertical line
  if (wt.critical_alpha) {
    const cx = xOf(wt.critical_alpha);
    ctx.strokeStyle = "rgba(245, 158, 11, 0.7)"; ctx.setLineDash([6, 4]);
    ctx.beginPath(); ctx.moveTo(cx, pad.t); ctx.lineTo(cx, H - pad.b); ctx.stroke(); ctx.setLineDash([]);
    ctx.fillStyle = "#fbbf24"; ctx.font = "bold 11px monospace";
    ctx.fillText(`α_c = ${wt.critical_alpha} /s (Jam Onset)`, Math.min(W - 220, cx + 6), pad.t + 14);
  }

  const drawCurve = (pts, yFn, color, width, dashed = false) => {
    ctx.strokeStyle = color; ctx.lineWidth = width; ctx.setLineDash(dashed ? [5, 4] : []);
    ctx.beginPath();
    pts.forEach((p, idx) => {
      const x = xOf(p.alpha), y = yFn(p);
      if (idx === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    });
    ctx.stroke(); ctx.setLineDash([]);
    pts.forEach((p) => {
      ctx.fillStyle = color;
      ctx.beginPath(); ctx.arc(xOf(p.alpha), yFn(p), 3.6, 0, Math.PI * 2); ctx.fill();
    });
  };

  drawCurve(wt.points, (p) => yFlux(p.protein_flux_per_min), "#2dd4bf", 2.6);
  drawCurve(wt.points, (p) => yNgd(p.ngd_risk_score), "#f43f5e", 2.2, true);
  if (opt) {
    drawCurve(opt.points, (p) => yFlux(p.protein_flux_per_min), "#38bdf8", 2.6);
  }

  ctx.fillStyle = "#2dd4bf"; ctx.font = "11px monospace";
  ctx.fillText(`${maxFlux.toFixed(0)} /min`, 6, pad.t + 8);
  ctx.fillText("Flux J", 8, pad.t + 24);
  ctx.fillStyle = "#f43f5e";
  ctx.fillText("100 NGD", W - 52, pad.t + 8);
  ctx.fillStyle = "#8ea3c4";
  ctx.fillText("Initiation Rate α (events / second)", W / 2 - 110, H - 10);
}

/* ---------- PRINTABLE LAB REPORT EXPORT ---------- */
function exportPrintableReport() {
  if (!state.simulation || !state.analysis) return;
  const sim = state.simulation;
  const an = state.analysis;
  const f5 = sim.folding_5prime;
  const win = window.open("", "_blank");
  if (!win) return;
  win.document.write(`<!DOCTYPE html><html><head><title>RiboFlow v2.0 Lab Report — ${state.currentName}</title>
  <style>
    body { font-family: Segoe UI, system-ui, sans-serif; color: #0f172a; max-width: 820px; margin: 28px auto; padding: 0 20px; line-height: 1.5; }
    h1 { border-bottom: 2px solid #0d9488; padding-bottom: 8px; }
    .grid { display: grid; grid-template-columns: repeat(3, 1fr); gap: 10px; margin: 16px 0; }
    .box { border: 1px solid #cbd5e1; border-radius: 8px; padding: 10px 12px; }
    .box span { display: block; font-size: 12px; color: #475569; }
    .box strong { font-size: 18px; }
    pre { background: #f8fafc; border: 1px solid #cbd5e1; padding: 12px; border-radius: 8px; overflow-wrap: anywhere; white-space: pre-wrap; font-size: 12px; }
  </style></head><body>
  <h1>RiboFlow v2.0 — Biophysical Translation &amp; mRNA Folding Report</h1>
  <p><strong>Transcript:</strong> ${state.currentName} | <strong>Host:</strong> ${sim.organism} | <strong>Length:</strong> ${an.codon_count} codons (${an.sequence.length} nt)</p>
  <div class="grid">
    <div class="box"><span>Protein Flux (J)</span><strong>${sim.protein_flux_per_min} /min</strong></div>
    <div class="box"><span>Effective α (after 5' ΔG)</span><strong>${sim.effective_initiation_rate} /s</strong></div>
    <div class="box"><span>5' Hairpin MFE (ΔG_5')</span><strong>${f5.mfe_kcal_mol} kcal/mol (${(f5.initiation_accessibility * 100).toFixed(0)}%)</strong></div>
    <div class="box"><span>Disome Collision %</span><strong>${(sim.collision_fraction * 100).toFixed(2)}% (${sim.collision_events} events)</strong></div>
    <div class="box"><span>NGD / RQC Risk</span><strong>${sim.ngd_risk_score} / 100</strong></div>
    <div class="box"><span>Peak Exit-Tunnel Charge</span><strong>+${sim.peak_tunnel_charge} K/R</strong></div>
  </div>
  <h3>5' mRNA Secondary Structure (${f5.window_nt} nt)</h3>
  <pre>RNA : ${f5.rna_sequence}\nFold: ${f5.dot_bracket}\nVerdict: ${f5.severity}</pre>
  <h3>Coding Sequence (CDS)</h3>
  <pre>>${state.currentName}\n${an.sequence}</pre>
  <p style="font-size:12px;color:#64748b;">Generated by RiboFlow v2.0 (https://emiresti.github.io/riboflow/)</p>
  <script>window.print();</script></body></html>`);
  win.document.close();
}

init();
