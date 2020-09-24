#![allow(unused_imports)]
use algebra_core::{
    curves::{models::SWModelParameters, AffineCurve, PairingEngine, ProjectiveCurve},
    fields::{Field, FpParameters, PrimeField, SquareRootField},
    test_rng, CanonicalSerialize, One, Zero,
};

use crate::bls12_377::{
    g1, g2, Bls12_377, Fq, Fq12, Fq2, Fr, G1Affine, G1Projective, G2Affine, G2Projective,
    Parameters,
};

use core::ops::{AddAssign, MulAssign, Neg};
use rand::Rng;

use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::prelude::*;

const NUM_TESTS: usize = 100;
const PREFIX: &str = "bls12377";
const FE_SIZE: usize = 48;
const WORD_SIZE: usize = 64;

#[derive(Serialize, Deserialize)]
struct VectorSuccess {
    input: String,
    expected: String,
}

fn write_vectors(vectors: Vec<VectorSuccess>, name: &str) {
    let serialized: String = serde_json::to_string(&vectors).unwrap();
    let mut file = File::create(PREFIX.to_string() + name + ".json").expect("must create the file");
    file.write(serialized.as_bytes())
        .expect("must write vectors");
}

#[test]
fn generate_test_vectors() {
    gen_g1_add_vectors();
    gen_g1_mul_vectors();
    gen_g1_multiexp_vectors();
    gen_g2_add_vectors();
    gen_g2_mul_vectors();
    gen_g2_multiexp_vectors();
    gen_pairing_vectors();
}

fn encode_g1(p: G1Projective) -> Vec<u8> {
    let fe_size: usize = 48;
    let word_size: usize = 64;
    let pad_zeros: Vec<u8> = vec![0u8; word_size - fe_size];
    let mut bytes: Vec<u8> = vec![];

    let mut buf_x = vec![];
    let p_affine = p.into_affine();

    p_affine
        .x
        .serialize(&mut buf_x)
        .expect("x coordinate must be serialized");
    bytes.extend(pad_zeros.clone());
    bytes.extend(buf_x.iter().rev());

    let mut buf_y = vec![];

    p_affine
        .y
        .serialize(&mut buf_y)
        .expect("y coordinate must be serialized");
    bytes.extend(pad_zeros.clone());
    bytes.extend(buf_y.iter().rev());

    bytes
}

fn encode_g2(p: G2Projective) -> Vec<u8> {
    let mut bytes: Vec<u8> = vec![];
    let pad_zeros: Vec<u8> = vec![0u8; WORD_SIZE - FE_SIZE];

    let mut buf = vec![];
    let p_affine = p.into_affine();

    p_affine
        .x
        .c0
        .serialize(&mut buf)
        .expect("c0 of x coordinate must be serialized");
    bytes.extend(pad_zeros.clone());
    bytes.extend(buf.iter().rev());
    buf.clear();

    p_affine
        .x
        .c1
        .serialize(&mut buf)
        .expect("c1 of x coordinate must be serialized");
    bytes.extend(pad_zeros.clone());
    bytes.extend(buf.iter().rev());
    buf.clear();

    p_affine
        .y
        .c0
        .serialize(&mut buf)
        .expect("c0 of y coordinate must be serialized");
    bytes.extend(pad_zeros.clone());
    bytes.extend(buf.iter().rev());
    buf.clear();

    p_affine
        .y
        .c1
        .serialize(&mut buf)
        .expect("c1 of y coordinate must be serialized");
    bytes.extend(pad_zeros.clone());
    bytes.extend(buf.iter().rev());

    bytes
}

fn encode_fr(p: Fr) -> Vec<u8> {
    let mut bytes = vec![];
    let mut buf = vec![];
    p.serialize(&mut buf).expect("scalar must be serialized");
    bytes.extend(buf.iter().rev());

    bytes
}

fn gen_g1_add_vectors() {
    let mut rng = test_rng();
    let mut vectors: Vec<VectorSuccess> = vec![];
    for _ in 0..NUM_TESTS {
        let mut input_bytes: Vec<u8> = vec![];
        let mut a: G1Projective = rng.gen();
        let b: G1Projective = rng.gen();
        let a_bytes: Vec<u8> = encode_g1(a);
        let b_bytes: Vec<u8> = encode_g1(b);
        input_bytes.extend(a_bytes);
        input_bytes.extend(b_bytes);
        let input: String = hex::encode(input_bytes.clone());

        a.add_assign(b);
        let result_bytes: Vec<u8> = encode_g1(a);
        let result: String = hex::encode(result_bytes);
        let vector = VectorSuccess {
            input,
            expected: result,
        };
        vectors.push(vector);
    }
    write_vectors(vectors, "_g1_add");
}

fn gen_g1_mul_vectors() {
    let mut rng = test_rng();
    let mut vectors: Vec<VectorSuccess> = vec![];
    for _ in 0..NUM_TESTS {
        let mut input_bytes: Vec<u8> = vec![];

        let mut a: G1Projective = rng.gen();
        let e: Fr = rng.gen();
        let a_bytes = encode_g1(a);
        let e_bytes = encode_fr(e);

        input_bytes.extend(a_bytes);
        input_bytes.extend(e_bytes);
        let input: String = hex::encode(input_bytes.clone());

        a.mul_assign(e);
        let result_bytes: Vec<u8> = encode_g1(a);
        let result: String = hex::encode(result_bytes);
        let vector = VectorSuccess {
            input,
            expected: result,
        };
        vectors.push(vector);
    }
    write_vectors(vectors, "_g1_mul");
}

fn gen_g1_multiexp_vectors() {
    let mut rng = test_rng();
    let mut vectors: Vec<VectorSuccess> = vec![];
    let mul_pair_size: usize = NUM_TESTS;
    for i in 1..mul_pair_size + 1 {
        let mut input_bytes: Vec<u8> = vec![];
        let mut acc = G1Projective::zero();
        for _ in 0..i {
            let mut a: G1Projective = rng.gen();
            let e: Fr = rng.gen();
            let a_bytes = encode_g1(a);
            let e_bytes = encode_fr(e);

            input_bytes.extend(a_bytes);
            input_bytes.extend(e_bytes);

            a.mul_assign(e);
            acc.add_assign(a);
        }
        let input: String = hex::encode(input_bytes.clone());

        let result_bytes: Vec<u8> = encode_g1(acc);
        let result: String = hex::encode(result_bytes);
        let vector = VectorSuccess {
            input,
            expected: result,
        };
        vectors.push(vector);
    }
    write_vectors(vectors, "_g1_multi_exp");
}

fn gen_g2_add_vectors() {
    let mut rng = test_rng();
    let mut vectors: Vec<VectorSuccess> = vec![];
    for _ in 0..NUM_TESTS {
        let mut input_bytes: Vec<u8> = vec![];
        let mut a: G2Projective = rng.gen();
        let b: G2Projective = rng.gen();
        let a_bytes: Vec<u8> = encode_g2(a);
        let b_bytes: Vec<u8> = encode_g2(b);
        input_bytes.extend(a_bytes);
        input_bytes.extend(b_bytes);
        let input: String = hex::encode(input_bytes.clone());

        a.add_assign(b);
        let result_bytes: Vec<u8> = encode_g2(a);
        let result: String = hex::encode(result_bytes);
        let vector = VectorSuccess {
            input,
            expected: result,
        };
        vectors.push(vector);
    }
    write_vectors(vectors, "_g2_add");
}

fn gen_g2_mul_vectors() {
    let mut rng = test_rng();
    let mut vectors: Vec<VectorSuccess> = vec![];
    for _ in 0..NUM_TESTS {
        let mut input_bytes: Vec<u8> = vec![];

        let mut a: G2Projective = rng.gen();
        let e: Fr = rng.gen();
        let a_bytes = encode_g2(a);
        let e_bytes = encode_fr(e);

        input_bytes.extend(a_bytes);
        input_bytes.extend(e_bytes);
        let input: String = hex::encode(input_bytes.clone());

        a.mul_assign(e);
        let result_bytes: Vec<u8> = encode_g2(a);
        let result: String = hex::encode(result_bytes);
        let vector = VectorSuccess {
            input,
            expected: result,
        };
        vectors.push(vector);
    }
    write_vectors(vectors, "_g2_mul");
}

fn gen_g2_multiexp_vectors() {
    let mut rng = test_rng();
    let mut vectors: Vec<VectorSuccess> = vec![];
    let mul_pair_size: usize = NUM_TESTS;
    for i in 1..mul_pair_size + 1 {
        let mut input_bytes: Vec<u8> = vec![];
        let mut acc = G2Projective::zero();
        for _ in 0..i {
            let mut a: G2Projective = rng.gen();
            let e: Fr = rng.gen();
            let a_bytes = encode_g2(a);
            let e_bytes = encode_fr(e);

            input_bytes.extend(a_bytes);
            input_bytes.extend(e_bytes);

            a.mul_assign(e);
            acc.add_assign(a);
        }
        let input: String = hex::encode(input_bytes.clone());

        let result_bytes: Vec<u8> = encode_g2(acc);
        let result: String = hex::encode(result_bytes);
        let vector = VectorSuccess {
            input,
            expected: result,
        };
        vectors.push(vector);
    }
    write_vectors(vectors, "_g2_multi_exp");
}

fn gen_pairing_vectors() {
    let mut rng = test_rng();
    let mut vectors: Vec<VectorSuccess> = vec![];
    let mut positive_result_bytes: Vec<u8> = vec![0u8; 32];
    positive_result_bytes[31] = 1u8;
    let negative_result_bytes: Vec<u8> = vec![0u8; 32];
    let g1_inf_encoded: Vec<u8> = vec![0u8; 128];
    let g2_inf_encoded: Vec<u8> = vec![0u8; 256];

    let g1 = G1Projective::prime_subgroup_generator();
    let g2 = G2Projective::prime_subgroup_generator();

    // expect true
    {
        // a. single pair
        {
            let mut input_bytes: Vec<u8> = vec![];

            let mut bytes_a1 = g1_inf_encoded.clone();
            let mut bytes_a2 = encode_g2(g2.clone());
            input_bytes.extend(bytes_a1);
            input_bytes.extend(bytes_a2);

            let input: String = hex::encode(input_bytes.clone());

            let vector = VectorSuccess {
                input,
                expected: hex::encode(positive_result_bytes.clone()),
            };
            vectors.push(vector);

            input_bytes.clear();
            bytes_a1 = encode_g1(g1.clone());
            bytes_a2 = g2_inf_encoded.to_vec().clone();
            input_bytes.extend(bytes_a1);
            input_bytes.extend(bytes_a2);

            let input: String = hex::encode(input_bytes.clone());

            let vector = VectorSuccess {
                input,
                expected: hex::encode(positive_result_bytes.clone()),
            };
            vectors.push(vector);
        }

        // b. multiple pair
        {
            for i in 0..NUM_TESTS {
                let mut acc: Fr = Fr::zero();
                let pair_size: usize = i + 2;
                let mut input_bytes: Vec<u8> = vec![];
                // n-1 pairs
                for _ in 0..pair_size - 1 {
                    let mut e1: Fr = rng.gen();
                    let e2: Fr = rng.gen();
                    let a1 = g1.mul(e1);
                    let a2 = g2.mul(e2);
                    let bytes_a1 = encode_g1(a1);
                    let bytes_a2 = encode_g2(a2);
                    input_bytes.extend(bytes_a1);
                    input_bytes.extend(bytes_a2);
                    // println!("e1\n{}", e1);
                    // println!("e2\n{}", e2);
                    // println!("acc\n{}", acc);
                    e1.mul_assign(e2);
                    acc.add_assign(e1);
                }
                // println!("acc\n{}", acc);
                // last pair
                let a1 = g1.mul(acc.neg());
                // println!("nacc\n{}", acc.neg());
                let a2 = g2.clone();
                let bytes_a1 = encode_g1(a1);
                let bytes_a2 = encode_g2(a2);
                input_bytes.extend(bytes_a1);
                input_bytes.extend(bytes_a2);

                let input: String = hex::encode(input_bytes.clone());
                let result: String = hex::encode(positive_result_bytes.clone());

                let vector = VectorSuccess {
                    input,
                    expected: result,
                };
                vectors.push(vector);
            }
        }
    }

    // expect false
    {
        for i in 0..NUM_TESTS {
            let pair_size: usize = i + 1;
            let mut input_bytes: Vec<u8> = vec![];
            for _ in 0..pair_size {
                let e1: Fr = rng.gen();
                let e2: Fr = rng.gen();
                let a1 = g1.mul(e1);
                let a2 = g2.mul(e2);
                let bytes_a1 = encode_g1(a1);
                let bytes_a2 = encode_g2(a2);
                input_bytes.extend(bytes_a1);
                input_bytes.extend(bytes_a2);
            }

            let input: String = hex::encode(input_bytes.clone());
            let result: String = hex::encode(negative_result_bytes.clone());

            let vector = VectorSuccess {
                input,
                expected: result,
            };
            vectors.push(vector);
        }
    }

    write_vectors(vectors, "_pairing");
}
