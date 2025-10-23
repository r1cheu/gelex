#include "gelex/data/genotype_loader.h"

#include <gelex/barkeep.h>

#include "gelex/data/genotype_pipe.h"
#include "gelex/logger.h"

namespace gelex
{

namespace bk = barkeep;

auto GenotypeLoader::load_from_bed(
    const std::filesystem::path& bed_path,
    std::shared_ptr<SampleManager> sample_manager,
    bool dominant,
    size_t chunk_size) -> std::expected<GenotypeMatrix, Error>
{
    return load_with_processor<NonStandardizingProcessor>(
        bed_path, std::move(sample_manager), dominant, chunk_size);
}

}  // namespace gelex
