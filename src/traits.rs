// Copyright (c) 2021 Toposware, Inc.
//
// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE or
// http://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or http://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.

use core::fmt::Debug;

// TODO: To support different kind of fields, we'd need them to
// implement some Field Trait which we can refer to here.

/// Defines output type of a cryptographic hash function.
pub trait Digest: Debug + Default + Copy + Clone + Eq + PartialEq + Send + Sync {
    /// Returns this digest serialized into an array of bytes.
    fn as_bytes(&self) -> [u8; 32];
}

/// Trait for implementing a cryptographic hash function.
pub trait Hasher {
    /// Specifies a digest type returned by this hasher.
    type Digest: Digest;

    /// Returns a hash of the provided sequence of bytes.
    fn hash(bytes: &[u8]) -> Self::Digest;

    /// Returns a hash of two digests.
    /// This method is intended for use in construction of Merkle trees.
    fn merge(values: &[Self::Digest; 2]) -> Self::Digest;

    /// Returns hash(`seed` || `value`).
    /// This method is intended for use in random coin sampling contexts.
    fn merge_with_int(seed: Self::Digest, value: u64) -> Self::Digest;
}
