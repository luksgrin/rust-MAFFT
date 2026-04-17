/// Dump BranchWeights for a 6-sequence topology to validate against C.

use mafft_tree::{Topology, BranchWeights, musclesupg, ClusterMethod, DistanceMatrix, JoinStep};

#[test]
fn dump_6seq_branch_weights() {
    let nseq = 6;
    let mut dm = DistanceMatrix::new(nseq);
    for i in 0..nseq {
        for j in (i + 1)..nseq {
            dm.set(i, j, (j - i) as f64 * 0.1);
        }
    }
    let topo = musclesupg(&dm, ClusterMethod::default());

    eprintln!("Topology ({nseq} seqs):");
    for (k, step) in topo.steps.iter().enumerate() {
        eprintln!("  step {k}: left={:?} right={:?} ll={:.6} rl={:.6}",
            step.left, step.right, step.left_length, step.right_length);
    }

    let bw = BranchWeights::new(&topo);

    for k in 0..topo.steps.len() {
        for side in 0..2 {
            let w = bw.weights_for_branch(&topo, k, side);
            let w_str: Vec<String> = w.iter().map(|v| format!("{:.8}", v)).collect();
            eprintln!("  w[step={k},side={side}] = [{}]", w_str.join(", "));
        }
    }

    // Sanity checks
    for k in 0..topo.steps.len() {
        for side in 0..2 {
            let w = bw.weights_for_branch(&topo, k, side);
            assert!(w.iter().all(|&v| v > 0.0 && v.is_finite()),
                "step={k} side={side}: bad weights {:?}", w);
        }
    }
}

#[test]
fn symmetric_4seq_weights_match_expected() {
    // ((0,1):0.1,(2,3):0.1):0.2 — perfectly symmetric
    let mut topo = Topology::new(4);
    topo.steps.push(JoinStep {
        left: vec![0], right: vec![1],
        left_length: 0.1, right_length: 0.1,
    });
    topo.steps.push(JoinStep {
        left: vec![2], right: vec![3],
        left_length: 0.1, right_length: 0.1,
    });
    topo.steps.push(JoinStep {
        left: vec![0, 1], right: vec![2, 3],
        left_length: 0.2, right_length: 0.2,
    });

    let bw = BranchWeights::new(&topo);

    // Step 2 (root), side 0: split {0,1} vs {2,3}
    // In a perfectly symmetric tree, ALL sequences should get the SAME weight.
    let w = bw.weights_for_branch(&topo, 2, 0);
    eprintln!("symmetric root split w = {:?}", w);
    let first = w[0];
    for (i, &v) in w.iter().enumerate() {
        assert!((v - first).abs() < 1e-10,
            "seq {i}: weight {v} != {first} in symmetric tree");
    }

    // Step 0 (leaf pair), side 0: split {0} vs {1,2,3}
    // seq 0 should be 1.0 (the split side, not modified)
    // seq 1 should differ from seqs 2,3 (sibling vs distant)
    let w = bw.weights_for_branch(&topo, 0, 0);
    eprintln!("leaf split w = {:?}", w);
    assert!((w[0] - 1.0).abs() < 1e-10, "split seq should be 1.0");
    assert!((w[2] - w[3]).abs() < 1e-10, "symmetric pair should match");
    // seq 1 is sibling (through node 0), seqs 2,3 are through node 0 → node 1
    // seq 1 should have higher weight than seqs 2,3
    assert!(w[1] > w[2], "sibling should have higher weight than distant: {} vs {}", w[1], w[2]);
}
