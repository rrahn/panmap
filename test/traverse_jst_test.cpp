// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <gtest/gtest.h>
#include <loadjst.hpp>
#include <libspm/matcher/horspool_matcher.hpp>
#include <libjst/sequence_tree/volatile_tree.hpp>
#include <libjst/sequence_tree/labelled_tree.hpp>
#include <libjst/sequence_tree/coloured_tree.hpp>
#include <libjst/sequence_tree/trim_tree.hpp>
#include <libjst/sequence_tree/prune_tree.hpp>
#include <libjst/sequence_tree/merge_tree.hpp>
#include <libjst/sequence_tree/seekable_tree.hpp>
#include <libjst/traversal/tree_traverser_base.hpp>
#include <libjst/sequence_tree/left_extend_tree.hpp>


TEST(test, traverse){
    using namespace spm;
    auto rcsstorepath = std::filesystem::path{std::string{DATADIR}}.concat("test.jst");
    rcs_store_t store = loadjst(rcsstorepath);
    auto query = "ACGC"_dna5;
    spm::horspool_matcher matcher{query};

    if (spm::window_size(matcher) == 0)
        return;

    auto search_tree = libjst::make_volatile(store) | libjst::labelled()
                                | libjst::coloured()
                                | libjst::trim(spm::window_size(matcher) - 1)
                                | libjst::prune()
                                | libjst::left_extend(spm::window_size(matcher) - 1)
                                | libjst::merge() // make big nodes
                                | libjst::seek();

    libjst::tree_traverser_base oblivious_path{search_tree};
    for (auto it = oblivious_path.begin(); it != oblivious_path.end(); ++it) {
        auto && cargo = *it;
        matcher(cargo.sequence(), [&] ([[maybe_unused]] auto && label_finder) {
            std::cout<<"Found hit. Yayy." << cargo.position() <<"\n";
            // callback(query, match_position{.tree_position{cargo.position()},
                                        //    .label_offset{std::ranges::ssize(cargo.sequence()) - seqan2::endPosition(label_finder)}});
        });
    }

}

