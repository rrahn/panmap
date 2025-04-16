// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <gtest/gtest.h>
#include <loadjst.hpp>
#include <libspm/seqan/alphabet.hpp>

TEST(test, create_rcs) {
    using namespace spm;

    rcs_store_t rcs_store{"ACGTGACTGACTGACTAGCTACG"_dna5, 5};
   using variant_t = std::ranges::range_value_t<cms_t>;
   auto dom = rcs_store.variants().coverage_domain();
   rcs_store.add (variant_t{libjst::breakpoint{3,4},"C"_dna5, coverage_t{{0,2},dom}});
   EXPECT_EQ(rcs_store.size(), 5);
}

TEST(test, loadrcsstore){
    auto rcsstorepath = std::filesystem::path{std::string{DATADIR}}.concat("test.jst");
    rcs_store_t store = loadjst(rcsstorepath);
    EXPECT_EQ(store.size(), 5);
}
