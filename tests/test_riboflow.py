"""Comprehensive unit test suite for RiboFlow v2.0."""

from __future__ import annotations

import unittest

from riboflow import (
    ORGANISMS,
    SAMPLES,
    analyze_sequence,
    apply_synonymous_mutations,
    compare_harmonization,
    compute_exit_tunnel_profile,
    design_synbio_sequence,
    find_restriction_hits,
    fold_5prime_window,
    get_organism,
    harmonize_codons,
    naive_optimize_codons,
    parse_fasta_or_raw,
    scan_critical_silent_mutations,
    simulate_ribosome_traffic,
    sweep_phase_diagram,
    translate_codons,
    unfold_5prime_synonymously,
    validate_cds,
)
from riboflow.web.server import get_metadata_payload, handle_api_post


class TestOrganismsAndSamples(unittest.TestCase):
    def test_all_organisms_have_64_codons(self) -> None:
        for key, org in ORGANISMS.items():
            self.assertEqual(len(org.codon_freq_per_1000), 64, msg=key)
            rates = org.elongation_rates()
            self.assertEqual(len(rates), 64)
            self.assertTrue(all(r > 0 for r in rates.values()))

    def test_all_curated_samples_are_valid_cds(self) -> None:
        self.assertGreaterEqual(len(SAMPLES), 6)
        for sample_id, sample in SAMPLES.items():
            codons = validate_cds(sample.sequence)
            self.assertGreaterEqual(len(codons), 50, msg=sample_id)
            self.assertEqual(codons[0], "ATG", msg=sample_id)
            self.assertIn(codons[-1], {"TAA", "TAG", "TGA"}, msg=sample_id)


class TestSequenceAnalyzer(unittest.TestCase):
    def test_fasta_and_rna_normalization(self) -> None:
        raw = ">TestGene\nAUGGCUGCUAAAGCUGCUGCUAAAGCUGCUGCUAAAGCUGCUGCUAAAUAA\n"
        header, dna = parse_fasta_or_raw(raw)
        self.assertEqual(header, "TestGene")
        self.assertNotIn("U", dna)
        codons = validate_cds(dna)
        self.assertEqual(codons[0], "ATG")

    def test_invalid_cds_rejected(self) -> None:
        with self.assertRaises(ValueError):
            validate_cds("ATGCG")
        with self.assertRaises(ValueError):
            validate_cds("GCT" * 20)
        with self.assertRaises(ValueError):
            validate_cds("ATG" + "GCT" * 10 + "TAA" + "GCT" * 10)

    def test_analysis_metrics_bounds(self) -> None:
        ecoli = get_organism("ecoli")
        res = analyze_sequence(SAMPLES["human_ins"].sequence, ecoli)
        self.assertGreater(res.cai, 0.0)
        self.assertLessEqual(res.cai, 1.0)
        self.assertGreater(res.gc_percent, 40.0)
        self.assertEqual(len(res.codons), res.codon_count)


class TestFivePrimeFolding(unittest.TestCase):
    def test_gc_hairpin_scores_stronger_than_open_au(self) -> None:
        tight = SAMPLES["polybasic_arrest"].sequence
        open_seq = "ATG" + "AAAAATAAACAA" * 6 + "TAA"
        rep_tight = fold_5prime_window(tight, window_nt=45)
        rep_open = fold_5prime_window(open_seq, window_nt=45)

        self.assertLess(rep_tight.mfe_kcal_mol, rep_open.mfe_kcal_mol)
        self.assertLess(rep_tight.initiation_accessibility, rep_open.initiation_accessibility)
        self.assertGreater(len(rep_tight.arcs), 0)

    def test_synonymous_5prime_unfolding_preserves_protein_and_improves_mfe(self) -> None:
        seq = SAMPLES["polybasic_arrest"].sequence
        res = unfold_5prime_synonymously(seq)
        self.assertTrue(res.protein_preserved)
        self.assertGreaterEqual(res.after.mfe_kcal_mol, res.before.mfe_kcal_mol)
        self.assertGreaterEqual(
            res.after.initiation_accessibility, res.before.initiation_accessibility
        )


class TestTASEPSimulatorAndPhysics(unittest.TestCase):
    def test_steric_exclusion_never_violated(self) -> None:
        ecoli = get_organism("ecoli")
        sim = simulate_ribosome_traffic(
            SAMPLES["jam_benchmark"].sequence,
            ecoli,
            initiation_rate=1.6,
            duration_seconds=80.0,
            capture_frames=80,
            seed=7,
        )
        self.assertGreater(sim.completed_proteins, 0)
        self.assertGreater(sim.collision_fraction, 0.0)
        for frame in sim.frames:
            positions = sorted(r.pos for r in frame.ribosomes)
            for a, b in zip(positions, positions[1:]):
                self.assertGreaterEqual(b - a, ecoli.ribosome_footprint)

    def test_exit_tunnel_polybasic_drag_detected(self) -> None:
        codons = validate_cds(SAMPLES["polybasic_arrest"].sequence)
        charges, drags = compute_exit_tunnel_profile(codons)
        self.assertGreater(max(charges), 6.0)
        self.assertLess(min(drags), 0.6)

        ecoli = get_organism("ecoli")
        sim_drag = simulate_ribosome_traffic(
            SAMPLES["polybasic_arrest"].sequence,
            ecoli,
            initiation_rate=1.4,
            couple_5prime_folding=False,
            enable_tunnel_drag=True,
            duration_seconds=80.0,
            capture_frames=0,
            seed=42,
        )
        sim_nodrag = simulate_ribosome_traffic(
            SAMPLES["polybasic_arrest"].sequence,
            ecoli,
            initiation_rate=1.4,
            couple_5prime_folding=False,
            enable_tunnel_drag=False,
            duration_seconds=80.0,
            capture_frames=0,
            seed=42,
        )
        self.assertLess(min(sim_drag.codon_rates), min(sim_nodrag.codon_rates))

    def test_phase_diagram_sweep(self) -> None:
        ecoli = get_organism("ecoli")
        rep = sweep_phase_diagram(
            SAMPLES["jam_benchmark"].sequence,
            ecoli,
            alpha_values=[0.1, 0.5, 1.2, 2.0],
            duration_seconds=45.0,
            seed=42,
        )
        self.assertEqual(len(rep.points), 4)
        self.assertGreater(rep.max_flux_per_min, 0.0)
        self.assertIsNotNone(rep.critical_alpha)

    def test_deterministic_with_seed(self) -> None:
        human = get_organism("human")
        s1 = simulate_ribosome_traffic(SAMPLES["cftr_nbd1"].sequence, human, seed=123)
        s2 = simulate_ribosome_traffic(SAMPLES["cftr_nbd1"].sequence, human, seed=123)
        self.assertEqual(s1.completed_proteins, s2.completed_proteins)
        self.assertEqual(s1.collision_events, s2.collision_events)
        self.assertEqual(s1.codon_occupancy, s2.codon_occupancy)

    def test_bottleneck_causes_more_collisions_than_fast_optimized(self) -> None:
        ecoli = get_organism("ecoli")
        jam_seq = SAMPLES["jam_benchmark"].sequence
        fast_seq = "".join(naive_optimize_codons(validate_cds(jam_seq), ecoli))

        jam_sim = simulate_ribosome_traffic(jam_seq, ecoli, initiation_rate=1.4, seed=42)
        fast_sim = simulate_ribosome_traffic(fast_seq, ecoli, initiation_rate=1.4, seed=42)

        self.assertGreater(jam_sim.ngd_risk_score, fast_sim.ngd_risk_score)
        self.assertGreater(jam_sim.collision_fraction, fast_sim.collision_fraction)


class TestSynBioAndHarmonizer(unittest.TestCase):
    def test_harmonization_preserves_exact_protein_and_improves_rhythm(self) -> None:
        human = get_organism("human")
        ecoli = get_organism("ecoli")
        seq = SAMPLES["human_ins"].sequence

        report = compare_harmonization(seq, source=human, host=ecoli, duration_seconds=60.0)
        self.assertTrue(report.protein_preserved)
        self.assertEqual(
            translate_codons(validate_cds(report.harmonized.sequence)),
            report.protein_sequence,
        )
        self.assertGreater(
            report.harmonized.rhythm_correlation,
            report.naive_optimized.rhythm_correlation,
        )

    def test_synbio_scrubs_restriction_sites_and_melts_hairpin(self) -> None:
        human = get_organism("human")
        ecoli = get_organism("ecoli")
        seq = SAMPLES["polybasic_arrest"].sequence
        self.assertGreaterEqual(len(find_restriction_hits(seq)), 3)

        rep = design_synbio_sequence(
            seq,
            source=human,
            host=ecoli,
            preserve_pauses=True,
            unfold_5prime=True,
            scrub_restriction_sites=True,
            avoid_homopolymers=True,
        )
        self.assertTrue(rep.protein_preserved)
        self.assertEqual(len(rep.after_qc.restriction_hits), 0)
        self.assertGreater(
            rep.after_qc.folding_5prime.mfe_kcal_mol,
            rep.before_qc.folding_5prime.mfe_kcal_mol,
        )

    def test_mrna_vaccine_uridine_depletion_lowers_u_content(self) -> None:
        human = get_organism("human")
        seq = SAMPLES["sars_cov2_rbd"].sequence
        rep = design_synbio_sequence(
            seq,
            source=human,
            host=human,
            preserve_pauses=True,
            unfold_5prime=True,
            uridine_deplete=True,
        )
        self.assertTrue(rep.protein_preserved)
        self.assertLess(rep.after_qc.uridine_percent, rep.before_qc.uridine_percent - 5.0)
        self.assertLess(rep.after_qc.uu_dinucleotides, rep.before_qc.uu_dinucleotides)

    def test_harmonize_codons_handles_all_samples(self) -> None:
        for sample in SAMPLES.values():
            src = get_organism(sample.organism)
            host = get_organism(sample.recommended_host)
            wt = validate_cds(sample.sequence)
            harm = harmonize_codons(wt, src, host)
            self.assertEqual(translate_codons(wt), translate_codons(harm))


class TestSilentMutationLab(unittest.TestCase):
    def test_synonymous_mutation_enforced(self) -> None:
        ecoli = get_organism("ecoli")
        seq = SAMPLES["human_ins"].sequence
        codons = validate_cds(seq)
        target_idx = next(i for i in range(5, len(codons) - 1) if codons[i] == "CTG")
        comp = apply_synonymous_mutations(seq, ecoli, edits={target_idx: "CTA"}, duration_seconds=60.0)
        self.assertTrue(comp.protein_identical)
        self.assertEqual(len(comp.edits), 1)
        with self.assertRaises(ValueError):
            apply_synonymous_mutations(seq, ecoli, edits={target_idx: "AAA"})

    def test_hotspot_scanner_returns_valid_silent_swaps(self) -> None:
        ecoli = get_organism("ecoli")
        hotspots = scan_critical_silent_mutations(
            SAMPLES["cftr_nbd1"].sequence,
            ecoli,
            initiation_rate=1.3,
            top_k=4,
            seed=42,
        )
        self.assertGreater(len(hotspots), 0)
        for h in hotspots:
            self.assertEqual(h.nucleotide_changes, 1)
            self.assertGreater(h.rate_drop_percent, 0.0)


class TestWebAPIRouterAndCLI(unittest.TestCase):
    def test_all_v2_api_endpoints(self) -> None:
        meta = get_metadata_payload()
        self.assertEqual(len(meta["organisms"]), 4)
        self.assertGreaterEqual(len(meta["samples"]), 6)

        st_sim, b_sim = handle_api_post(
            "/api/simulate",
            {"sequence": SAMPLES["polybasic_arrest"].sequence, "organism": "ecoli", "duration": 40},
        )
        self.assertEqual(st_sim, 200)
        self.assertIn("folding_5prime", b_sim["simulation"])
        self.assertIn("tunnel_charge_profile", b_sim["simulation"])

        st_fold, b_fold = handle_api_post("/api/fold", {"sequence": SAMPLES["polybasic_arrest"].sequence})
        self.assertEqual(st_fold, 200)
        self.assertIn("arcs", b_fold)

        st_unf, b_unf = handle_api_post("/api/unfold", {"sequence": SAMPLES["polybasic_arrest"].sequence})
        self.assertEqual(st_unf, 200)
        self.assertTrue(b_unf["protein_preserved"])

        st_syn, b_syn = handle_api_post(
            "/api/synbio",
            {
                "sequence": SAMPLES["sars_cov2_rbd"].sequence,
                "source": "human",
                "host": "human",
                "uridine_deplete": True,
            },
        )
        self.assertEqual(st_syn, 200)
        self.assertTrue(b_syn["synbio"]["protein_preserved"])

        st_ph, b_ph = handle_api_post(
            "/api/phase",
            {"sequence": SAMPLES["jam_benchmark"].sequence, "organism": "ecoli"},
        )
        self.assertEqual(st_ph, 200)
        self.assertIn("wildtype", b_ph)

        st_err, _ = handle_api_post("/api/simulate", {"sequence": "NOT_DNA_123", "organism": "ecoli"})
        self.assertEqual(st_err, 400)

    def test_cli_v2_subcommands(self) -> None:
        from riboflow.cli import main

        self.assertEqual(main(["analyze", "--sample", "human_ins", "-o", "ecoli"]), 0)
        self.assertEqual(main(["fold", "--sample", "polybasic_arrest", "--melt"]), 0)
        self.assertEqual(main(["simulate", "--sample", "polybasic_arrest", "-o", "ecoli", "-d", "40"]), 0)
        self.assertEqual(main(["synbio", "--sample", "sars_cov2_rbd", "-s", "human", "-t", "human", "--vaccine-u-deplete"]), 0)
        self.assertEqual(main(["phase", "--sample", "jam_benchmark", "-o", "ecoli"]), 0)


if __name__ == "__main__":
    unittest.main()
