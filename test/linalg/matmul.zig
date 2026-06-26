const std = @import("std");

const zsl = @import("zsl");

const tzsl = @import("../zsl.zig");

const combinations = blk: {
    @setEvalBranchQuota(292_200);
    break :blk [_][2]type{
        // matgensta_matgensta
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },

        // matgensta_matsymsta
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },

        // matgensta_matsymden
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major) },

        // matgensta_matsymspa
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major) },

        // matgensta_mathersta
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },

        // matgensta_matherden
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major) },

        // matgensta_matherspa
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major) },

        // matgensta_mattrista
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },

        // matgensta_matdiasta
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), tzsl.matrix.diagonal.Static(zsl.cf64) },

        // matgensta_matpersta
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },

        // matgensta_matperspa
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), zsl.matrix.permutation.Sparse(zsl.cf64, .forward) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), zsl.matrix.permutation.Sparse(zsl.cf64, .backward) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), zsl.matrix.permutation.Sparse(zsl.cf64, .forward) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), zsl.matrix.permutation.Sparse(zsl.cf64, .backward) },

        // matgensta_vecsta
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), tzsl.vector.Static(zsl.cf64) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), tzsl.vector.Static(zsl.cf64) },

        // matgensta_vecden
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), zsl.vector.Dense(zsl.cf64) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), zsl.vector.Dense(zsl.cf64) },

        // matgensta_vecspa
        .{ tzsl.matrix.general.Static(zsl.cf64, .row_major), zsl.vector.Sparse(zsl.cf64) },
        .{ tzsl.matrix.general.Static(zsl.cf64, .col_major), zsl.vector.Sparse(zsl.cf64) },

        // matsymsta_matgensta
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },

        // matsymsta_matsymsta
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },

        // matsymsta_matsymden
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major) },

        // matsymsta_matsymspa
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major) },

        // matsymsta_mathersta
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },

        // matsymsta_matherden
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major) },

        // matsymsta_matherspa
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major) },

        // matsymsta_mattrista
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },

        // matsymsta_matdiasta
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.diagonal.Static(zsl.cf64) },

        // matsymsta_matpersta
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },

        // matsymsta_matperspa
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), zsl.matrix.permutation.Sparse(zsl.cf64, .forward) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), zsl.matrix.permutation.Sparse(zsl.cf64, .backward) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), zsl.matrix.permutation.Sparse(zsl.cf64, .forward) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), zsl.matrix.permutation.Sparse(zsl.cf64, .backward) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), zsl.matrix.permutation.Sparse(zsl.cf64, .forward) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), zsl.matrix.permutation.Sparse(zsl.cf64, .backward) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), zsl.matrix.permutation.Sparse(zsl.cf64, .forward) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), zsl.matrix.permutation.Sparse(zsl.cf64, .backward) },

        // matsymsta_vecsta
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), tzsl.vector.Static(zsl.cf64) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), tzsl.vector.Static(zsl.cf64) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), tzsl.vector.Static(zsl.cf64) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), tzsl.vector.Static(zsl.cf64) },

        // matsymsta_vecden
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), zsl.vector.Dense(zsl.cf64) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), zsl.vector.Dense(zsl.cf64) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), zsl.vector.Dense(zsl.cf64) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), zsl.vector.Dense(zsl.cf64) },

        // matsymsta_vecspa
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major), zsl.vector.Sparse(zsl.cf64) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major), zsl.vector.Sparse(zsl.cf64) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major), zsl.vector.Sparse(zsl.cf64) },
        .{ tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major), zsl.vector.Sparse(zsl.cf64) },

        // matsymden_matgensta
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },

        // matsymden_matsymsta
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },

        // matsymden_mathersta
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },

        // matsymden_mattrista
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },

        // matsymden_matdiasta
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.diagonal.Static(zsl.cf64) },

        // matsymden_matpersta
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },

        // matsymden_vecsta
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major), tzsl.vector.Static(zsl.cf64) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major), tzsl.vector.Static(zsl.cf64) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major), tzsl.vector.Static(zsl.cf64) },
        .{ zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major), tzsl.vector.Static(zsl.cf64) },

        // matsymspa_matgensta
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },

        // matsymspa_matsymsta
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },

        // matsymspa_mathersta
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },

        // matsymspa_mattrista
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },

        // matsymspa_matdiasta
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.diagonal.Static(zsl.cf64) },

        // matsymspa_matpersta
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },

        // matsymspa_vecsta
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major), tzsl.vector.Static(zsl.cf64) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major), tzsl.vector.Static(zsl.cf64) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major), tzsl.vector.Static(zsl.cf64) },
        .{ zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major), tzsl.vector.Static(zsl.cf64) },

        // mathersta_matgensta
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },

        // mathersta_matsymsta
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },

        // mathersta_matsymden
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major) },

        // mathersta_matsymspa
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major) },

        // mathersta_mathersta
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },

        // mathersta_matherden
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major) },

        // mathersta_matherspa
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major) },

        // mathersta_mattrista
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },

        // mathersta_matdiasta
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.diagonal.Static(zsl.cf64) },

        // mathersta_matpersta
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },

        // mathersta_matperspa
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), zsl.matrix.permutation.Sparse(zsl.cf64, .forward) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), zsl.matrix.permutation.Sparse(zsl.cf64, .backward) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), zsl.matrix.permutation.Sparse(zsl.cf64, .forward) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), zsl.matrix.permutation.Sparse(zsl.cf64, .backward) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), zsl.matrix.permutation.Sparse(zsl.cf64, .forward) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), zsl.matrix.permutation.Sparse(zsl.cf64, .backward) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), zsl.matrix.permutation.Sparse(zsl.cf64, .forward) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), zsl.matrix.permutation.Sparse(zsl.cf64, .backward) },

        // mathersta_vecsta
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), tzsl.vector.Static(zsl.cf64) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), tzsl.vector.Static(zsl.cf64) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), tzsl.vector.Static(zsl.cf64) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), tzsl.vector.Static(zsl.cf64) },

        // mathersta_vecden
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), zsl.vector.Dense(zsl.cf64) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), zsl.vector.Dense(zsl.cf64) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), zsl.vector.Dense(zsl.cf64) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), zsl.vector.Dense(zsl.cf64) },

        // mathersta_vecspa
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major), zsl.vector.Sparse(zsl.cf64) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major), zsl.vector.Sparse(zsl.cf64) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major), zsl.vector.Sparse(zsl.cf64) },
        .{ tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major), zsl.vector.Sparse(zsl.cf64) },

        // matherden_matgensta
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },

        // matherden_matsymsta
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },

        // matherden_mathersta
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },

        // matherden_mattrista
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },

        // matherden_matdiasta
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.diagonal.Static(zsl.cf64) },

        // matherden_matpersta
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },

        // matherden_vecsta
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major), tzsl.vector.Static(zsl.cf64) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major), tzsl.vector.Static(zsl.cf64) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major), tzsl.vector.Static(zsl.cf64) },
        .{ zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major), tzsl.vector.Static(zsl.cf64) },

        // matherspa_matgensta
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },

        // matherspa_matsymsta
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },

        // matherspa_mathersta
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },

        // matherspa_mattrista
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },

        // matherspa_matdiasta
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.diagonal.Static(zsl.cf64) },

        // matherspa_matpersta
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },

        // matherspa_vecsta
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major), tzsl.vector.Static(zsl.cf64) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major), tzsl.vector.Static(zsl.cf64) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major), tzsl.vector.Static(zsl.cf64) },
        .{ zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major), tzsl.vector.Static(zsl.cf64) },

        // mattrista_matgensta
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), tzsl.matrix.general.Static(zsl.cf64, .col_major) },

        // mattrista_matsymsta
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },

        // mattrista_matsymden
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major) },

        // mattrista_matsymspa
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major) },

        // mattrista_mathersta
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },

        // mattrista_matherden
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major) },

        // mattrista_matherspa
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major) },

        // mattrista_mattrista
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },

        // mattrista_matdiasta
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), tzsl.matrix.diagonal.Static(zsl.cf64) },

        // mattrista_matpersta
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },

        // mattrista_matperspa
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), zsl.matrix.permutation.Sparse(zsl.cf64, .forward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), zsl.matrix.permutation.Sparse(zsl.cf64, .backward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), zsl.matrix.permutation.Sparse(zsl.cf64, .forward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), zsl.matrix.permutation.Sparse(zsl.cf64, .backward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), zsl.matrix.permutation.Sparse(zsl.cf64, .forward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), zsl.matrix.permutation.Sparse(zsl.cf64, .backward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), zsl.matrix.permutation.Sparse(zsl.cf64, .forward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), zsl.matrix.permutation.Sparse(zsl.cf64, .backward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), zsl.matrix.permutation.Sparse(zsl.cf64, .forward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), zsl.matrix.permutation.Sparse(zsl.cf64, .backward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), zsl.matrix.permutation.Sparse(zsl.cf64, .forward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), zsl.matrix.permutation.Sparse(zsl.cf64, .backward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), zsl.matrix.permutation.Sparse(zsl.cf64, .forward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), zsl.matrix.permutation.Sparse(zsl.cf64, .backward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), zsl.matrix.permutation.Sparse(zsl.cf64, .forward) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), zsl.matrix.permutation.Sparse(zsl.cf64, .backward) },

        // mattrista_vecsta
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), tzsl.vector.Static(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), tzsl.vector.Static(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), tzsl.vector.Static(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), tzsl.vector.Static(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), tzsl.vector.Static(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), tzsl.vector.Static(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), tzsl.vector.Static(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), tzsl.vector.Static(zsl.cf64) },

        // mattrista_vecden
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), zsl.vector.Dense(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), zsl.vector.Dense(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), zsl.vector.Dense(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), zsl.vector.Dense(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), zsl.vector.Dense(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), zsl.vector.Dense(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), zsl.vector.Dense(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), zsl.vector.Dense(zsl.cf64) },

        // mattrista_vecspa
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major), zsl.vector.Sparse(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major), zsl.vector.Sparse(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major), zsl.vector.Sparse(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major), zsl.vector.Sparse(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major), zsl.vector.Sparse(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major), zsl.vector.Sparse(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major), zsl.vector.Sparse(zsl.cf64) },
        .{ tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major), zsl.vector.Sparse(zsl.cf64) },

        // matdiasta_matgensta
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), tzsl.matrix.general.Static(zsl.cf64, .col_major) },

        // matdiasta_matsymsta
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },

        // matdiasta_matsymden
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major) },

        // matdiasta_matsymspa
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major) },

        // matdiasta_mathersta
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },

        // matdiasta_matherden
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major) },

        // matdiasta_matherspa
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major) },

        // matdiasta_mattrista
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },

        // matdiasta_matdiasta
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), tzsl.matrix.diagonal.Static(zsl.cf64) },

        // matdiasta_matpersta
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },

        // matdiasta_matperspa
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), zsl.matrix.permutation.Sparse(zsl.cf64, .forward) },
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), zsl.matrix.permutation.Sparse(zsl.cf64, .backward) },

        // matdiasta_vecsta
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), tzsl.vector.Static(zsl.cf64) },

        // matdiasta_vecden
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), zsl.vector.Dense(zsl.cf64) },

        // matdiasta_vecspa
        .{ tzsl.matrix.diagonal.Static(zsl.cf64), zsl.vector.Sparse(zsl.cf64) },

        // matpersta_matgensta
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), tzsl.matrix.general.Static(zsl.cf64, .col_major) },

        // matpersta_matsymsta
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },

        // matpersta_matsymden
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major) },

        // matpersta_matsymspa
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major) },

        // matpersta_mathersta
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },

        // matpersta_matherden
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major) },

        // matpersta_matherspa
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major) },

        // matpersta_mattrista
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },

        // matpersta_matdiasta
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), tzsl.matrix.diagonal.Static(zsl.cf64) },

        // matpersta_matpersta
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },

        // matpersta_matperspa
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), zsl.matrix.permutation.Sparse(zsl.cf64, .forward) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), zsl.matrix.permutation.Sparse(zsl.cf64, .backward) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), zsl.matrix.permutation.Sparse(zsl.cf64, .forward) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), zsl.matrix.permutation.Sparse(zsl.cf64, .backward) },

        // matpersta_vecsta
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), tzsl.vector.Static(zsl.cf64) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), tzsl.vector.Static(zsl.cf64) },

        // matpersta_vecden
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), zsl.vector.Dense(zsl.cf64) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), zsl.vector.Dense(zsl.cf64) },

        // matpersta_vecspa
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .forward), zsl.vector.Sparse(zsl.cf64) },
        .{ tzsl.matrix.permutation.Static(zsl.cf64, .backward), zsl.vector.Sparse(zsl.cf64) },

        // matperspa_matgensta
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .forward), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .forward), tzsl.matrix.general.Static(zsl.cf64, .col_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .backward), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .backward), tzsl.matrix.general.Static(zsl.cf64, .col_major) },

        // matperspa_matsymsta
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .forward), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .forward), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .forward), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .forward), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .backward), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .backward), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .backward), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .backward), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },

        // matperspa_mathersta
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .forward), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .forward), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .forward), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .forward), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .backward), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .backward), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .backward), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .backward), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },

        // matperspa_mattrista
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .forward), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .forward), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .forward), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .forward), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .forward), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .forward), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .forward), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .forward), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .backward), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .backward), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .backward), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .backward), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .backward), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .backward), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .backward), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .backward), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },

        // matperspa_matdiasta
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .forward), tzsl.matrix.diagonal.Static(zsl.cf64) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .backward), tzsl.matrix.diagonal.Static(zsl.cf64) },

        // matperspa_matpersta
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .forward), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .forward), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .backward), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .backward), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },

        // matperspa_vecsta
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .forward), tzsl.vector.Static(zsl.cf64) },
        .{ zsl.matrix.permutation.Sparse(zsl.cf64, .backward), tzsl.vector.Static(zsl.cf64) },

        // vecsta_matgensta
        .{ tzsl.vector.Static(zsl.cf64), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ tzsl.vector.Static(zsl.cf64), tzsl.matrix.general.Static(zsl.cf64, .col_major) },

        // vecsta_matsymsta
        .{ tzsl.vector.Static(zsl.cf64), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.vector.Static(zsl.cf64), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.vector.Static(zsl.cf64), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.vector.Static(zsl.cf64), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },

        // vecsta_matsymden
        .{ tzsl.vector.Static(zsl.cf64), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.vector.Static(zsl.cf64), zsl.matrix.symmetric.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.vector.Static(zsl.cf64), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.vector.Static(zsl.cf64), zsl.matrix.symmetric.Dense(zsl.cf64, .lower, .col_major) },

        // vecsta_matsymspa
        .{ tzsl.vector.Static(zsl.cf64), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.vector.Static(zsl.cf64), zsl.matrix.symmetric.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.vector.Static(zsl.cf64), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.vector.Static(zsl.cf64), zsl.matrix.symmetric.Sparse(zsl.cf64, .lower, .col_major) },

        // vecsta_mathersta
        .{ tzsl.vector.Static(zsl.cf64), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ tzsl.vector.Static(zsl.cf64), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ tzsl.vector.Static(zsl.cf64), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ tzsl.vector.Static(zsl.cf64), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },

        // vecsta_matherden
        .{ tzsl.vector.Static(zsl.cf64), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .row_major) },
        .{ tzsl.vector.Static(zsl.cf64), zsl.matrix.hermitian.Dense(zsl.cf64, .upper, .col_major) },
        .{ tzsl.vector.Static(zsl.cf64), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .row_major) },
        .{ tzsl.vector.Static(zsl.cf64), zsl.matrix.hermitian.Dense(zsl.cf64, .lower, .col_major) },

        // vecsta_matherspa
        .{ tzsl.vector.Static(zsl.cf64), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .row_major) },
        .{ tzsl.vector.Static(zsl.cf64), zsl.matrix.hermitian.Sparse(zsl.cf64, .upper, .col_major) },
        .{ tzsl.vector.Static(zsl.cf64), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .row_major) },
        .{ tzsl.vector.Static(zsl.cf64), zsl.matrix.hermitian.Sparse(zsl.cf64, .lower, .col_major) },

        // vecsta_mattrista
        .{ tzsl.vector.Static(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ tzsl.vector.Static(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ tzsl.vector.Static(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ tzsl.vector.Static(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ tzsl.vector.Static(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ tzsl.vector.Static(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ tzsl.vector.Static(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ tzsl.vector.Static(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },

        // vecsta_matdiasta
        .{ tzsl.vector.Static(zsl.cf64), tzsl.matrix.diagonal.Static(zsl.cf64) },

        // vecsta_matpersta
        .{ tzsl.vector.Static(zsl.cf64), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ tzsl.vector.Static(zsl.cf64), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },

        // vecsta_matperspa
        .{ tzsl.vector.Static(zsl.cf64), zsl.matrix.permutation.Sparse(zsl.cf64, .forward) },
        .{ tzsl.vector.Static(zsl.cf64), zsl.matrix.permutation.Sparse(zsl.cf64, .backward) },

        // vecden_matgensta
        .{ zsl.vector.Dense(zsl.cf64), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ zsl.vector.Dense(zsl.cf64), tzsl.matrix.general.Static(zsl.cf64, .col_major) },

        // vecden_matsymsta
        .{ zsl.vector.Dense(zsl.cf64), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.vector.Dense(zsl.cf64), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.vector.Dense(zsl.cf64), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.vector.Dense(zsl.cf64), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },

        // vecden_mathersta
        .{ zsl.vector.Dense(zsl.cf64), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.vector.Dense(zsl.cf64), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.vector.Dense(zsl.cf64), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.vector.Dense(zsl.cf64), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },

        // vecden_mattrista
        .{ zsl.vector.Dense(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ zsl.vector.Dense(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ zsl.vector.Dense(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ zsl.vector.Dense(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ zsl.vector.Dense(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ zsl.vector.Dense(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ zsl.vector.Dense(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ zsl.vector.Dense(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },

        // vecden_matdiasta
        .{ zsl.vector.Dense(zsl.cf64), tzsl.matrix.diagonal.Static(zsl.cf64) },

        // vecden_matpersta
        .{ zsl.vector.Dense(zsl.cf64), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ zsl.vector.Dense(zsl.cf64), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },

        // vecspa_matgensta
        .{ zsl.vector.Sparse(zsl.cf64), tzsl.matrix.general.Static(zsl.cf64, .row_major) },
        .{ zsl.vector.Sparse(zsl.cf64), tzsl.matrix.general.Static(zsl.cf64, .col_major) },

        // vecspa_matsymsta
        .{ zsl.vector.Sparse(zsl.cf64), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.vector.Sparse(zsl.cf64), tzsl.matrix.symmetric.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.vector.Sparse(zsl.cf64), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.vector.Sparse(zsl.cf64), tzsl.matrix.symmetric.Static(zsl.cf64, .lower, .col_major) },

        // vecspa_mathersta
        .{ zsl.vector.Sparse(zsl.cf64), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .row_major) },
        .{ zsl.vector.Sparse(zsl.cf64), tzsl.matrix.hermitian.Static(zsl.cf64, .upper, .col_major) },
        .{ zsl.vector.Sparse(zsl.cf64), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .row_major) },
        .{ zsl.vector.Sparse(zsl.cf64), tzsl.matrix.hermitian.Static(zsl.cf64, .lower, .col_major) },

        // vecspa_mattrista
        .{ zsl.vector.Sparse(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .row_major) },
        .{ zsl.vector.Sparse(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .non_unit, .col_major) },
        .{ zsl.vector.Sparse(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .row_major) },
        .{ zsl.vector.Sparse(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .upper, .unit, .col_major) },
        .{ zsl.vector.Sparse(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .row_major) },
        .{ zsl.vector.Sparse(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .non_unit, .col_major) },
        .{ zsl.vector.Sparse(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .row_major) },
        .{ zsl.vector.Sparse(zsl.cf64), tzsl.matrix.triangular.Static(zsl.cf64, .lower, .unit, .col_major) },

        // vecspa_matdiasta
        .{ zsl.vector.Sparse(zsl.cf64), tzsl.matrix.diagonal.Static(zsl.cf64) },

        // vecspa_matpersta
        .{ zsl.vector.Sparse(zsl.cf64), tzsl.matrix.permutation.Static(zsl.cf64, .forward) },
        .{ zsl.vector.Sparse(zsl.cf64), tzsl.matrix.permutation.Static(zsl.cf64, .backward) },
    };
};

const mkn_limits = [_][3]usize{
    .{ 7, 7, 7 },
    .{ 7, 16, 7 },
    .{ 7, 16, 16 },
    .{ 16, 7, 7 },
    .{ 16, 16, 16 },
    .{ 33, 33, 16 },
    .{ 33, 7, 16 },
    .{ 33, 16, 7 },
    .{ 7, 33, 16 },
};

test "zsl.matrix.matmul" {
    @setEvalBranchQuota(1_000_000);

    const allocator = std.testing.allocator;

    var prng = std.Random.DefaultPrng.init(@bitCast(std.Io.Clock.real.now(std.testing.io).toSeconds()));
    const rand = prng.random();

    inline for (combinations) |fake_combination| {
        inline for (mkn_limits) |mkn| {
            comptime var combination = fake_combination;
            inline for (0..combination.len) |fake_i| {
                const i = combination.len - 1 - fake_i;

                const m = if (comptime i == 0) mkn[0] else mkn[1];
                const n = if (comptime i == 0) mkn[1] else mkn[2];
                if (@hasDecl(combination[i], "instantiate")) {
                    if (@hasDecl(combination[i], "general"))
                        combination[i] = zsl.matrix.general.Static(m, n, combination[i].Numeric, combination[i].storage_layout)
                    else if (@hasDecl(combination[i], "symmetric"))
                        combination[i] = zsl.matrix.symmetric.Static(m, combination[i].Numeric, combination[i].storage_uplo, combination[i].storage_layout)
                    else if (@hasDecl(combination[i], "hermitian"))
                        combination[i] = zsl.matrix.hermitian.Static(m, combination[i].Numeric, combination[i].storage_uplo, combination[i].storage_layout)
                    else if (@hasDecl(combination[i], "triangular"))
                        combination[i] = zsl.matrix.triangular.Static(m, n, combination[i].Numeric, combination[i].storage_uplo, combination[i].storage_diag, combination[i].storage_layout)
                    else if (@hasDecl(combination[i], "diagonal"))
                        combination[i] = zsl.matrix.diagonal.Static(m, n, combination[i].Numeric)
                    else if (@hasDecl(combination[i], "permutation"))
                        combination[i] = zsl.matrix.permutation.Static(m, combination[i].Numeric, combination[i].storage_direction)
                    else
                        combination[i] = zsl.vector.Static(mkn[1], combination[i].Numeric);
                }
            }

            if (comptime !((zsl.meta.isSquareMatrix(combination[0]) and mkn[0] != mkn[1]) or
                (zsl.meta.isSquareMatrix(combination[1]) and mkn[1] != mkn[2])))
            {
                try executeTestBlock(allocator, rand, combination, mkn[0], mkn[1], mkn[2]);
            }
        }
    }
}

fn executeTestBlock(
    allocator: std.mem.Allocator,
    rand: std.Random,
    comptime combination: [2]type,
    comptime m: usize,
    comptime k: usize,
    comptime n: usize,
) !void {
    var B = if (comptime zsl.meta.isVector(combination[0]))
        try tzsl.vector.randomVector(
            combination[0],
            allocator,
            rand,
            k,
            1,
            k / 7,
        )
    else
        try tzsl.matrix.randomMatrix(
            combination[0],
            allocator,
            rand,
            m,
            k,
            zsl.int.min(m, k),
        );
    defer tzsl.deinit(allocator, &B);

    var C = if (comptime zsl.meta.isVector(combination[1]))
        try tzsl.vector.randomVector(
            combination[1],
            allocator,
            rand,
            k,
            1,
            k / 7,
        )
    else
        try tzsl.matrix.randomMatrix(
            combination[1],
            allocator,
            rand,
            k,
            n,
            zsl.int.min(k, n),
        );
    defer tzsl.deinit(allocator, &C);

    var D = try tzsl.linalg.correctMatmul(zsl.cf64, allocator, m, k, n, B, C);
    defer D.deinit(allocator);

    const A = zsl.linalg.matmul(B, C);

    if (comptime zsl.meta.isMatrix(@TypeOf(A))) {
        tzsl.matrix.areEql(A, D) catch |e| {
            std.debug.print("Failed on A: {s} = B: {s} * C: {s}, case m = {}, k = {}, n = {}\n", .{ @typeName(@TypeOf(A)), @typeName(combination[0]), @typeName(combination[1]), m, k, n });

            tzsl.matrix.printMatrix("A", A);
            tzsl.matrix.printMatrix("B", B);
            tzsl.matrix.printMatrix("C", C);
            tzsl.matrix.printMatrix("D", D);

            return e;
        };
    } else {
        tzsl.vector.areEql(A, D) catch |e| {
            std.debug.print("Failed on A: {s} = B: {s} * C: {s}, case m = {}, k = {}, n = {}\n", .{ @typeName(@TypeOf(A)), @typeName(combination[0]), @typeName(combination[1]), m, k, n });

            tzsl.vector.printVector("A", A);
            if (comptime zsl.meta.isMatrix(@TypeOf(B))) tzsl.matrix.printMatrix("B", B) else tzsl.vector.printVector("B", B);
            if (comptime zsl.meta.isMatrix(@TypeOf(C))) tzsl.matrix.printMatrix("C", C) else tzsl.vector.printVector("C", C);
            tzsl.vector.printVector("D", D);

            return e;
        };
    }
}
