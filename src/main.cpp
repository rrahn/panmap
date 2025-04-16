// // SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// // SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// // SPDX-License-Identifier: CC0-1.0

#include <sharg/all.hpp>
#include <fstream>
#include <span>
#include <ranges>
#include <seqan3/alignment/cigar_conversion/cigar_from_alignment.hpp>
#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sam_file/output.hpp>
#include <seqan3/io/sequence_file/input.hpp>

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

#include <configuration.hpp>
#include <loadjst.hpp>
#include <seqan3/search/search.hpp>

std::vector<libjst::seek_position> jst_search(const rcs_store_t& jst_data, const reference_t& read)
{

    spm::horspool_matcher matcher{read};

    auto search_tree = libjst::make_volatile(jst_data) | libjst::labelled()
                                | libjst::coloured()
                                | libjst::trim(spm::window_size(matcher) - 1)
                                | libjst::prune()
                                | libjst::left_extend(spm::window_size(matcher) - 1)
                                | libjst::merge() // make big nodes
                                | libjst::seek();

    std::vector<libjst::seek_position> hits;
    libjst::tree_traverser_base oblivious_path{search_tree};
    for (auto it = oblivious_path.begin(); it != oblivious_path.end(); ++it) {
        auto && cargo = *it;
        matcher(cargo.sequence(), [&] ([[maybe_unused]] auto && label_finder) {
            hits.push_back(cargo.position());
            std::cout<<"Found hit. Yayy." << cargo.position() <<"\n";
            // callback(query, match_position{.tree_position{cargo.position()},
                                        //    .label_offset{std::ranges::ssize(cargo.sequence()) - seqan2::endPosition(label_finder)}});
        });
    }

    return hits;
}

void map_reads(std::filesystem::path const & query_path,
    std::filesystem::path const & sam_path,
    const rcs_store_t & jst_data,
    uint8_t const errors)
{
    sequence_file_t query_file_in{query_path};
    seqan3::sam_file_output sam_out{sam_path, seqan3::fields<seqan3::field::seq,
                                              seqan3::field::id,
                                              seqan3::field::ref_id,
                                              seqan3::field::ref_offset,
                                              seqan3::field::cigar,
                                              seqan3::field::qual,
                                              seqan3::field::mapq>{}};



    // definiert die Ausrichtungskonfiguration mit einem globalen Alignment-Algorithmus, keine gap penalties
    seqan3::configuration const align_config =
        seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                         seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
                                         seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                         seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}} |
        seqan3::align_cfg::edit_scheme |
        seqan3::align_cfg::output_alignment{} |
        seqan3::align_cfg::output_begin_position{} |
        seqan3::align_cfg::output_score{};

    for (auto && record : query_file_in) // beginnt Schleife, die jeden Eintrag in der Query-Datei durchläuft
    {
        auto & query = record.sequence(); // referenziert die Sequenz des aktuellen Eintrags
        auto hits = jst_search(jst_data, query);

        for (auto hit_pos : hits)
        {
            // std::span text_view{std::data(jst_data.sequences[hit_pos]), query.size() +1};
            std::cout << "seek position: " << hit_pos << "\n";

            // for (auto&& alignment : seqan3::align_pairwise(std::tie(text_view,query), align_config))
            // {
            // auto cigar = seqan3::cigar_from_alignment(alignment.alignment()); // wandelt das Alignment in einen CIGAR-String um
            // size_t ref_offset = alignment.sequence1_begin_position() + 2; // berechnet Offset in der Referenzsequenz basierend auf der Ausrichtung
            // size_t map_qual = 60u + alignment.score(); // berechnet die Mapping-Qualität basierend auf dem Ausrichtungsscore
            //     //fängt einen neuen Eintrag zur SAM-Datei hinzu
            //     sam_out.emplace_back(query,
            //                          record.id(),
            //                          jst_data.ids[hit_pos],
            //                          ref_offset,
            //                          cigar,
            //                          record.base_qualities(),
            //                          map_qual);
            // }
        }
    }
}
// Funktion, die das Programm ausfühhrt
void run_program(std::filesystem::path const & reference_path,
                 std::filesystem::path const & query_path,
                 std::filesystem::path const & sam_path,
                 uint8_t const errors)
{
    //JST Ladefunktion???
 auto jst_data = loadjst(reference_path);

map_reads(query_path,
          sam_path,
          jst_data,
          errors);
}

// Funktion, die den Argument-Parser konfiguriert
void initialise_argument_parser(sharg::parser & parser, configuration & args)
{
    parser.info.author = "Julia Bodnar"; // setzt den Autor des Programms
    parser.info.short_description = "JST-based read mapper"; // setzt Beschreibung des Programms
    parser.info.version = "1.0.0"; // setzt Version des Programms

    parser.add_option(args.reference_path, // fügt eine Option für den Pfad zur Referenzdatei hinzu
                      sharg::config{.short_id = 'r', .long_id="reference", .description = "referenzgenome (JST-format)", .required = true, .validator = sharg::input_file_validator{{"jst"}}});

    parser.add_option(args.query_path,
                      sharg::config{.short_id = 'q',
                                    .long_id = "query",
                                    .description = "Input reads",
                                    .required = true, .validator = sharg::input_file_validator{{"fq","fastq"}}});

    parser.add_option(args.sam_path, // fügt eine Option für den Pfad zur SAM Ausgabedatei hinzu
                      sharg::config{.short_id = 'o',
                                    .long_id = "output",
                                    .description = "The output SAM file.",
                                    .validator = sharg::output_file_validator{{"sam"}}});

    parser.add_option(args.errors, // fügt eine Option für die maximale Anzahl von Feldern hinzu
                      sharg::config{.short_id = 'e',
                                    .long_id = "error",
                                    .description = "Maximum allowed errors.",
                                    .validator = sharg::arithmetic_range_validator{0, 4}});
}

int main(int argc, char const** argv)
{
    sharg::parser parser{"PanMap", argc, argv};
    configuration args{};

    initialise_argument_parser(parser, args);

    try {
        parser.parse();
    } catch (const sharg::parser_error& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    run_program(args.reference_path,
               args.query_path,
               args.sam_path,
               args.errors);

    return 0;
}
