#![allow(unused_imports)]
use algebra_core::{
    curves::{models::SWModelParameters, AffineCurve, PairingEngine, ProjectiveCurve},
    fields::{Field, FpParameters, PrimeField, SquareRootField},
    test_rng, CanonicalSerialize, One, Zero,
};

use crate::bw6_761::*;

use core::ops::{AddAssign, Mul, MulAssign, Neg};
use rand::Rng;

use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::prelude::*;

const NUM_TESTS: usize = 100;
const PREFIX: &str = "bw6";

const EVM_WORD_SIZE: usize = 32;
const FE_SIZE: usize = 96;
const FE_WORD_SIZE: usize = FE_SIZE;
const FR_SIZE: usize = 48;
const FR_WORD_SIZE: usize = 2 * EVM_WORD_SIZE;
const FR_ZERO_OFFSET: usize = FR_WORD_SIZE - FR_SIZE;
const G_WORD_SIZE: usize = 2 * FE_WORD_SIZE;

#[derive(Serialize, Deserialize)]
struct VectorSuccess {
    input: String,
    expected: String,
    name: String,
}

#[derive(Serialize, Deserialize)]
struct VectorFail {
    input: String,
    expected_error: String,
    name: String,
}

fn write_vectors(vectors: Vec<VectorSuccess>, name: &str) {
    let serialized: String = serde_json::to_string(&vectors).unwrap();
    let mut file = File::create(PREFIX.to_string() + name + ".json").expect("must create the file");
    file.write(serialized.as_bytes())
        .expect("must write vectors");
}

fn write_vectors_fail(vectors: Vec<VectorFail>, name: &str) {
    let serialized: String = serde_json::to_string(&vectors).unwrap();
    let mut file = File::create(PREFIX.to_string() + name + ".json").expect("must create the file");
    file.write(serialized.as_bytes())
        .expect("must write vectors");
}

fn encoded_fe_larger_than_modulus() -> Vec<u8> {
    hex::decode("0122e824fb83ce0ad187c94004faff3eb926186a81d14688528275ef8087be41707ba638e584e91903cebaff25b423048689c8ed12f9fd9071dcd3dc73ebff2e98a116c25667a8f8160cf8aeeaf0a437e6913e6870000082f49d00000000008f")
        .expect("must decode")
}

fn rand_g1_point_not_on_correct_subgroup() -> G1Projective {
    let mut rng = test_rng();

    loop {
        let x: Fq = rng.gen();
        let mut y: Fq = x.mul(x);
        y.mul_assign(x);
        y.add_assign(g1::Parameters::COEFF_B);
        if let Some(y) = y.sqrt() {
            let p = G1Affine::new(x, y, false);
            assert!(p.is_on_curve());
            assert!(!p.is_in_correct_subgroup_assuming_on_curve());
            return p.into_projective();
        }
    }
}

fn rand_g2_point_not_on_correct_subgroup() -> G2Projective {
    let mut rng = test_rng();

    loop {
        let x: Fq = rng.gen();
        let mut y: Fq = x.mul(x);
        y.mul_assign(x);
        y.add_assign(g2::Parameters::COEFF_B);
        if let Some(y) = y.sqrt() {
            let p = G2Affine::new(x, y, false);
            assert!(p.is_on_curve());
            assert!(!p.is_in_correct_subgroup_assuming_on_curve());
            return p.into_projective();
        }
    }
}

fn rand_g1_point_not_on_curve() -> G1Projective {
    let mut rng = test_rng();
    let x: Fq = rng.gen();
    let y: Fq = rng.gen();
    let p = G1Affine::new(x, y, false);
    assert!(!p.is_on_curve());
    p.into_projective()
}

fn rand_g2_point_not_on_curve() -> G2Projective {
    let mut rng = test_rng();
    let x: Fq = rng.gen();
    let y: Fq = rng.gen();
    let p = G2Affine::new(x, y, false);
    assert!(!p.is_on_curve());
    p.into_projective()
}

#[test]
fn generate_test_vectors() {
    // gen_g1_add_vectors();
    // gen_g1_mul_vectors();
    // gen_g1_multiexp_vectors();
    // gen_g2_add_vectors();
    // gen_g2_mul_vectors();
    // gen_g2_multiexp_vectors();
    // gen_pairing_vectors();
    gen_fail_g1_add_vectors();
    gen_fail_g1_mul_vectors();
    gen_fail_g1_multiexp_vectors();
    gen_fail_g2_add_vectors();
    gen_fail_g2_mul_vectors();
    gen_fail_g2_multiexp_vectors();
    gen_fail_pairing();
}

fn encode_g1(p: G1Projective) -> Vec<u8> {
    let mut bytes: Vec<u8> = vec![];
    let mut buf_x = vec![];
    let p_affine = p.into_affine();

    p_affine
        .x
        .serialize(&mut buf_x)
        .expect("x coordinate must be serialized");
    bytes.extend(buf_x.iter().rev());

    let mut buf_y = vec![];

    p_affine
        .y
        .serialize(&mut buf_y)
        .expect("y coordinate must be serialized");
    bytes.extend(buf_y.iter().rev());

    bytes
}

fn encode_g2(p: G2Projective) -> Vec<u8> {
    let mut bytes: Vec<u8> = vec![];

    let mut buf = vec![];
    let p_affine = p.into_affine();

    p_affine
        .x
        .serialize(&mut buf)
        .expect("x coordinate must be serialized");
    bytes.extend(buf.iter().rev());
    buf.clear();

    p_affine
        .y
        .serialize(&mut buf)
        .expect("y coordinate must be serialized");
    bytes.extend(buf.iter().rev());
    buf.clear();

    bytes
}

fn encode_fr(p: Fr) -> Vec<u8> {
    let mut bytes = vec![];
    let pad_zeros: Vec<u8> = vec![0u8; FR_WORD_SIZE - FR_SIZE];
    let mut buf = vec![];
    p.serialize(&mut buf).expect("scalar must be serialized");
    bytes.extend(pad_zeros.clone());
    bytes.extend(buf.iter().rev());

    bytes
}

fn gen_g1_add_vectors() {
    let mut rng = test_rng();
    let mut vectors: Vec<VectorSuccess> = vec![];
    for i in 0..NUM_TESTS {
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
            name: format!("{}_{}", "g1_add", i + 1),
        };

        vectors.push(vector);
    }
    write_vectors(vectors, "_g1_add");
}

fn gen_g1_mul_vectors() {
    let mut rng = test_rng();
    let mut vectors: Vec<VectorSuccess> = vec![];
    for i in 0..NUM_TESTS {
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
            name: format!("{}_{}", "g1_mul", i + 1),
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
            name: format!("{}_{}", "g1_multiexp", i + 1),
        };
        vectors.push(vector);
    }
    write_vectors(vectors, "_g1_multi_exp");
}

fn gen_g2_add_vectors() {
    let mut rng = test_rng();
    let mut vectors: Vec<VectorSuccess> = vec![];
    for i in 0..NUM_TESTS {
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
            name: format!("{}_{}", "g2_add", i + 1),
        };
        vectors.push(vector);
    }
    write_vectors(vectors, "_g2_add");
}

fn gen_g2_mul_vectors() {
    let mut rng = test_rng();
    let mut vectors: Vec<VectorSuccess> = vec![];
    for i in 0..NUM_TESTS {
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
            name: format!("{}_{}", "g2_mul", i + 1),
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
            name: format!("{}_{}", "g2_multiexp", i),
        };
        vectors.push(vector);
    }
    write_vectors(vectors, "_g2_multi_exp");
}

fn gen_pairing_vectors() {
    let mut rng = test_rng();
    let mut vectors: Vec<VectorSuccess> = vec![];
    let mut positive_result_bytes: Vec<u8> = vec![0u8; EVM_WORD_SIZE];
    positive_result_bytes[31] = 1u8;
    let negative_result_bytes: Vec<u8> = vec![0u8; EVM_WORD_SIZE];
    let g1_inf_encoded: Vec<u8> = vec![0u8; 2 * FE_SIZE];
    let g2_inf_encoded: Vec<u8> = vec![0u8; 2 * FE_SIZE];

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
                name: format!("{}", "g2_pairing_1"),
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
                name: format!("{}", "g2_pairing_2"),
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
                    e1.mul_assign(e2);
                    acc.add_assign(e1);
                }
                // last pair
                let a1 = g1.mul(acc.neg());
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
                    name: format!("{}_{}", "g2_pairing", i + 2),
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
                name: format!("{}_{}", "g2_pairing_0", NUM_TESTS + i + 2),
            };
            vectors.push(vector);
        }
    }

    write_vectors(vectors, "_pairing");
}

fn gen_fail_vectors(input_len: usize) -> Vec<VectorFail> {
    let mut vectors: Vec<VectorFail> = vec![];

    // invalid length: empty
    {
        let input: String = hex::encode(vec![]);
        let vector = VectorFail {
            input: String::from(""),
            expected_error: String::from("invalid input length"),
            name: format!("invalid_input_length_empty"),
        };
        vectors.push(vector);
    }

    // invalid length: short
    {
        let input: String = hex::encode(vec![0u8; input_len - 1]);
        let vector = VectorFail {
            input: String::from(""),
            expected_error: String::from("invalid input length"),
            name: format!("invalid_input_length_short"),
        };
        vectors.push(vector);
    }

    // invalid length: long
    {
        let input: String = hex::encode(vec![1u8; input_len + 1]);
        let vector = VectorFail {
            input,
            expected_error: String::from("invalid input length"),
            name: format!("invalid_input_length_large"),
        };
        vectors.push(vector);
    }

    // violate top zeros
    // {
    //     let input: String = hex::encode(vec![1u8; input_len]);
    //     let vector = VectorFail {
    //         input,
    //         expected_error: String::from("invalid field element top bytes"),
    //         name: format!("violate_top_zero_bytes"),
    //     };
    //     vectors.push(vector);
    // }

    vectors
}

fn gen_fail_g1_add_vectors() {
    let mut rng = test_rng();

    let input_len = 4 * FE_WORD_SIZE;

    let mut vectors: Vec<VectorFail> = gen_fail_vectors(input_len);

    // larger than modulus
    {
        let a: G1Projective = rng.gen();

        let mut input_bytes: Vec<u8> = vec![];
        let a_bytes = encode_g1(a);
        input_bytes.extend(a_bytes);
        input_bytes.extend(encoded_fe_larger_than_modulus());
        input_bytes.extend(vec![0u8; FE_SIZE]);

        let input: String = hex::encode(input_bytes.clone());
        let vector = VectorFail {
            input,
            expected_error: String::from("must be less than modulus"),
            name: format!("large_field_element"),
        };
        vectors.push(vector);
    }

    // not on curve
    {
        let a: G1Projective = rng.gen();
        let b: G1Projective = rand_g1_point_not_on_curve();

        let a_bytes = encode_g1(a);
        let e_bytes = encode_g1(b);

        let mut input_bytes: Vec<u8> = vec![];
        input_bytes.extend(a_bytes);
        input_bytes.extend(e_bytes);

        let input: String = hex::encode(input_bytes.clone());
        let vector = VectorFail {
            input,
            expected_error: String::from("point is not on curve"),
            name: format!("point_not_on_curve"),
        };
        vectors.push(vector);
    }
    write_vectors_fail(vectors, "_g1_add_fail");
}

fn gen_fail_g1_mul_vectors() {
    let input_len = 2 * FE_WORD_SIZE + FR_WORD_SIZE;
    // let pad_zeros: Vec<u8> = vec![0u8; FR_ZERO_OFFSET];
    let mut vectors: Vec<VectorFail> = gen_fail_vectors(input_len);

    // large modulus
    {
        let mut input_bytes: Vec<u8> = vec![];
        // x
        input_bytes.extend(encoded_fe_larger_than_modulus());
        // y
        input_bytes.extend(vec![0u8; FE_WORD_SIZE]);
        // e
        input_bytes.extend(vec![0u8; FR_WORD_SIZE]);

        let input: String = hex::encode(input_bytes.clone());
        let vector = VectorFail {
            input,
            expected_error: String::from("must be less than modulus"),
            name: format!("large_field_element"),
        };
        vectors.push(vector);
    }

    // not on curve
    {
        let a: G1Projective = rand_g1_point_not_on_curve();

        let a_bytes = encode_g1(a);

        let mut input_bytes: Vec<u8> = vec![];
        input_bytes.extend(a_bytes);
        input_bytes.extend(vec![0u8; FR_WORD_SIZE]);

        let input: String = hex::encode(input_bytes.clone());
        let vector = VectorFail {
            input,
            expected_error: String::from("point is not on curve"),
            name: format!("point_not_on_curve"),
        };
        vectors.push(vector);
    }

    // TODO: violate top zeros of fr
    write_vectors_fail(vectors, "_g1_mul_fail");
}

fn gen_fail_g1_multiexp_vectors() {
    let mut rng = test_rng();

    let input_len = 3 * (2 * FE_WORD_SIZE + FR_WORD_SIZE);
    let pad_zeros: Vec<u8> = vec![0u8; FR_ZERO_OFFSET];

    let mut vectors: Vec<VectorFail> = gen_fail_vectors(input_len);

    // large modulus
    {
        let a: G1Projective = rng.gen();
        let e1: Fr = rng.gen();
        let b: G1Projective = rng.gen();
        let e2: Fr = rng.gen();

        let mut input_bytes: Vec<u8> = vec![];

        let a_bytes = encode_g1(a);
        let e1_bytes = encode_fr(e1);
        input_bytes.extend(a_bytes);
        input_bytes.extend(e1_bytes);

        let b_bytes = encode_g1(b);
        let e2_bytes = encode_fr(e2);
        input_bytes.extend(b_bytes);
        input_bytes.extend(e2_bytes);

        input_bytes.extend(encoded_fe_larger_than_modulus());
        // y
        input_bytes.extend(vec![0u8; FE_WORD_SIZE]);
        // e
        input_bytes.extend(vec![0u8; FR_WORD_SIZE]);

        let input: String = hex::encode(input_bytes.clone());
        let vector = VectorFail {
            input,
            expected_error: String::from("must be less than modulus"),
            name: format!("large_field_element"),
        };
        vectors.push(vector);
    }

    // not on curve
    {
        let a: G1Projective = rng.gen();
        let e1: Fr = rng.gen();
        let b: G1Projective = rng.gen();
        let e2: Fr = rng.gen();
        let c = rand_g1_point_not_on_curve();
        let e3: Fr = rng.gen();

        let mut input_bytes: Vec<u8> = vec![];

        let a_bytes = encode_g1(a);
        let e1_bytes = encode_fr(e1);
        input_bytes.extend(a_bytes);
        input_bytes.extend(e1_bytes);

        let b_bytes = encode_g1(b);
        let e2_bytes = encode_fr(e2);
        input_bytes.extend(b_bytes);
        input_bytes.extend(e2_bytes);

        let c_bytes = encode_g1(c);
        let e3_bytes = encode_fr(e3);
        input_bytes.extend(c_bytes);
        input_bytes.extend(e3_bytes);

        let input: String = hex::encode(input_bytes.clone());
        let vector = VectorFail {
            input,
            expected_error: String::from("point is not on curve"),
            name: format!("point_not_on_curve"),
        };
        vectors.push(vector);
    }
    write_vectors_fail(vectors, "_g1_multiexp_fail");
}

fn gen_fail_g2_add_vectors() {
    let mut rng = test_rng();

    let input_len = 4 * FE_WORD_SIZE;

    let mut vectors: Vec<VectorFail> = gen_fail_vectors(input_len);

    // larger than modulus
    {
        let a: G2Projective = rng.gen();

        let mut input_bytes: Vec<u8> = vec![];
        let a_bytes = encode_g2(a);
        input_bytes.extend(a_bytes);
        input_bytes.extend(encoded_fe_larger_than_modulus());
        input_bytes.extend(vec![0u8; FE_WORD_SIZE]);

        let input: String = hex::encode(input_bytes.clone());
        let vector = VectorFail {
            input,
            expected_error: String::from("must be less than modulus"),
            name: format!("large_field_element"),
        };
        vectors.push(vector);
    }

    // not on curve
    {
        let a: G2Projective = rng.gen();
        let b: G2Projective = rand_g2_point_not_on_curve();

        let a_bytes = encode_g2(a);
        let e_bytes = encode_g2(b);

        let mut input_bytes: Vec<u8> = vec![];
        input_bytes.extend(a_bytes);
        input_bytes.extend(e_bytes);

        let input: String = hex::encode(input_bytes.clone());
        let vector = VectorFail {
            input,
            expected_error: String::from("point is not on curve"),
            name: format!("point_not_on_curve"),
        };
        vectors.push(vector);
    }
    write_vectors_fail(vectors, "_g2_add_fail");
}

fn gen_fail_g2_mul_vectors() {
    let input_len = 2 * FE_WORD_SIZE + FR_WORD_SIZE;
    // let pad_zeros: Vec<u8> = vec![0u8; FR_ZERO_OFFSET];
    let mut vectors: Vec<VectorFail> = gen_fail_vectors(input_len);

    // large modulus
    {
        let mut input_bytes: Vec<u8> = vec![];
        // x
        input_bytes.extend(encoded_fe_larger_than_modulus());
        // y
        input_bytes.extend(vec![0u8; FE_WORD_SIZE]);
        // e
        input_bytes.extend(vec![0u8; FR_WORD_SIZE]);

        let input: String = hex::encode(input_bytes.clone());
        let vector = VectorFail {
            input,
            expected_error: String::from("must be less than modulus"),
            name: format!("large_field_element"),
        };
        vectors.push(vector);
    }

    // not on curve
    {
        let a: G2Projective = rand_g2_point_not_on_curve();

        let a_bytes = encode_g2(a);

        let mut input_bytes: Vec<u8> = vec![];
        input_bytes.extend(a_bytes);
        input_bytes.extend(vec![0u8; FR_WORD_SIZE]);

        let input: String = hex::encode(input_bytes.clone());
        let vector = VectorFail {
            input,
            expected_error: String::from("point is not on curve"),
            name: format!("point_not_on_curve"),
        };
        vectors.push(vector);
    }

    // TODO: violate top zeros of fr
    write_vectors_fail(vectors, "_g2_mul_fail");
}

fn gen_fail_g2_multiexp_vectors() {
    let mut rng = test_rng();

    let input_len = 3 * (2 * FE_WORD_SIZE + FR_WORD_SIZE);
    let pad_zeros: Vec<u8> = vec![0u8; FR_ZERO_OFFSET];

    let mut vectors: Vec<VectorFail> = gen_fail_vectors(input_len);

    // large modulus
    {
        let a: G2Projective = rng.gen();
        let e1: Fr = rng.gen();
        let b: G2Projective = rng.gen();
        let e2: Fr = rng.gen();

        let mut input_bytes: Vec<u8> = vec![];

        let a_bytes = encode_g2(a);
        let e1_bytes = encode_fr(e1);
        input_bytes.extend(a_bytes);
        input_bytes.extend(e1_bytes);

        let b_bytes = encode_g2(b);
        let e2_bytes = encode_fr(e2);
        input_bytes.extend(b_bytes);
        input_bytes.extend(e2_bytes);

        input_bytes.extend(encoded_fe_larger_than_modulus());
        // y
        input_bytes.extend(vec![0u8; FE_WORD_SIZE]);
        // e
        input_bytes.extend(vec![0u8; FR_WORD_SIZE]);

        let input: String = hex::encode(input_bytes.clone());
        let vector = VectorFail {
            input,
            expected_error: String::from("must be less than modulus"),
            name: format!("large_field_element"),
        };
        vectors.push(vector);
    }

    // not on curve
    {
        let a: G2Projective = rng.gen();
        let e1: Fr = rng.gen();
        let b: G2Projective = rng.gen();
        let e2: Fr = rng.gen();
        let c = rand_g2_point_not_on_curve();
        let e3: Fr = rng.gen();

        let mut input_bytes: Vec<u8> = vec![];

        let a_bytes = encode_g2(a);
        let e1_bytes = encode_fr(e1);
        input_bytes.extend(a_bytes);
        input_bytes.extend(e1_bytes);

        let b_bytes = encode_g2(b);
        let e2_bytes = encode_fr(e2);
        input_bytes.extend(b_bytes);
        input_bytes.extend(e2_bytes);

        let c_bytes = encode_g2(c);
        let e3_bytes = encode_fr(e3);
        input_bytes.extend(c_bytes);
        input_bytes.extend(e3_bytes);

        let input: String = hex::encode(input_bytes.clone());
        let vector = VectorFail {
            input,
            expected_error: String::from("point is not on curve"),
            name: format!("point_not_on_curve"),
        };
        vectors.push(vector);
    }
    write_vectors_fail(vectors, "_g2_multiexp_fail");
}

fn gen_fail_pairing() {
    let mut rng = test_rng();
    let input_len = 3 * 4 * G_WORD_SIZE;
    let mut vectors: Vec<VectorFail> = gen_fail_vectors(input_len);

    // large modulus
    {
        let mut input_bytes: Vec<u8> = vec![];

        let a1: G1Projective = rng.gen();
        let a2: G2Projective = rng.gen();
        let a1_bytes = encode_g1(a1);
        let a2_bytes = encode_g2(a2);
        input_bytes.extend(a1_bytes);
        input_bytes.extend(a2_bytes);

        let b1: G1Projective = rng.gen();
        let b2: G2Projective = rng.gen();
        let b1_bytes = encode_g1(b1);
        let b2_bytes = encode_g2(b2);
        input_bytes.extend(b1_bytes);
        input_bytes.extend(b2_bytes);

        input_bytes.extend(encoded_fe_larger_than_modulus());
        input_bytes.extend(vec![0u8; 3 * FE_WORD_SIZE]);

        let input: String = hex::encode(input_bytes.clone());
        let vector = VectorFail {
            input,
            expected_error: String::from("must be less than modulus"),
            name: format!("large_field_element"),
        };
        vectors.push(vector);
    }

    // not on curve g1
    {
        let mut input_bytes: Vec<u8> = vec![];

        let a1: G1Projective = rng.gen();
        let a2: G2Projective = rng.gen();
        let a1_bytes = encode_g1(a1);
        let a2_bytes = encode_g2(a2);
        input_bytes.extend(a1_bytes);
        input_bytes.extend(a2_bytes);

        let b1: G1Projective = rng.gen();
        let b2: G2Projective = rng.gen();
        let b1_bytes = encode_g1(b1);
        let b2_bytes = encode_g2(b2);
        input_bytes.extend(b1_bytes);
        input_bytes.extend(b2_bytes);

        let c1: G1Projective = rand_g1_point_not_on_curve();
        let c2: G2Projective = rng.gen();
        let c1_bytes = encode_g1(c1);
        let c2_bytes = encode_g2(c2);
        input_bytes.extend(c1_bytes);
        input_bytes.extend(c2_bytes);

        let input: String = hex::encode(input_bytes.clone());
        let vector = VectorFail {
            input,
            expected_error: String::from("point is not on curve"),
            name: format!("point_not_on_curve_g1"),
        };
        vectors.push(vector);
    }

    // not on curve g2
    {
        let mut input_bytes: Vec<u8> = vec![];

        let a1: G1Projective = rng.gen();
        let a2: G2Projective = rng.gen();
        let a1_bytes = encode_g1(a1);
        let a2_bytes = encode_g2(a2);
        input_bytes.extend(a1_bytes);
        input_bytes.extend(a2_bytes);

        let b1: G1Projective = rng.gen();
        let b2: G2Projective = rng.gen();
        let b1_bytes = encode_g1(b1);
        let b2_bytes = encode_g2(b2);
        input_bytes.extend(b1_bytes);
        input_bytes.extend(b2_bytes);

        let c1: G1Projective = rng.gen();
        let c2: G2Projective = rand_g2_point_not_on_curve();
        let c1_bytes = encode_g1(c1);
        let c2_bytes = encode_g2(c2);
        input_bytes.extend(c1_bytes);
        input_bytes.extend(c2_bytes);

        let input: String = hex::encode(input_bytes.clone());
        let vector = VectorFail {
            input,
            expected_error: String::from("point is not on curve"),
            name: format!("point_not_on_curve_g2"),
        };
        vectors.push(vector);
    }

    // incorrect subgroup g1
    {
        let mut input_bytes: Vec<u8> = vec![];

        let a1: G1Projective = rng.gen();
        let a2: G2Projective = rng.gen();
        let a1_bytes = encode_g1(a1);
        let a2_bytes = encode_g2(a2);
        input_bytes.extend(a1_bytes);
        input_bytes.extend(a2_bytes);

        let b1: G1Projective = rng.gen();
        let b2: G2Projective = rng.gen();
        let b1_bytes = encode_g1(b1);
        let b2_bytes = encode_g2(b2);
        input_bytes.extend(b1_bytes);
        input_bytes.extend(b2_bytes);

        let c1: G1Projective = rand_g1_point_not_on_correct_subgroup();
        let c2: G2Projective = rng.gen();
        let c1_bytes = encode_g1(c1);
        let c2_bytes = encode_g2(c2);
        input_bytes.extend(c1_bytes);
        input_bytes.extend(c2_bytes);

        let input: String = hex::encode(input_bytes.clone());
        let vector = VectorFail {
            input,
            expected_error: String::from("g1 point is not on correct subgroup"),
            name: format!("incorrect_subgroup_g1"),
        };
        vectors.push(vector);
    }

    // incorrect subgroup g2
    {
        let mut input_bytes: Vec<u8> = vec![];

        let a1: G1Projective = rng.gen();
        let a2: G2Projective = rng.gen();
        let a1_bytes = encode_g1(a1);
        let a2_bytes = encode_g2(a2);
        input_bytes.extend(a1_bytes);
        input_bytes.extend(a2_bytes);

        let b1: G1Projective = rng.gen();
        let b2: G2Projective = rng.gen();
        let b1_bytes = encode_g1(b1);
        let b2_bytes = encode_g2(b2);
        input_bytes.extend(b1_bytes);
        input_bytes.extend(b2_bytes);

        let c1: G1Projective = rng.gen();
        let c2: G2Projective = rand_g2_point_not_on_correct_subgroup();
        let c1_bytes = encode_g1(c1);
        let c2_bytes = encode_g2(c2);
        input_bytes.extend(c1_bytes);
        input_bytes.extend(c2_bytes);

        let input: String = hex::encode(input_bytes.clone());
        let vector = VectorFail {
            input,
            expected_error: String::from("g2 point is not on correct subgroup"),
            name: format!("incorrect_subgroup_g2"),
        };
        vectors.push(vector);
    }

    write_vectors_fail(vectors, "_g2_pairing_fail");
}
