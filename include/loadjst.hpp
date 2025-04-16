// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <cerrno>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>

#include <globaltypes.hpp>

struct JST_Data {
    std::vector<std::string> ids;
    std::vector<seqan3::dna5_vector> sequences;
};

rcs_store_t loadjst(std::filesystem::path const & path)
{
    rcs_store_t rcsstore{};
    std::ifstream rcsstream{path, std::ios::binary};
    {
        cereal::BinaryInputArchive rcsarchive{rcsstream};
        rcsstore.load(rcsarchive);
    }
    return rcsstore;
}

