#include "bed_fixture.h"

#include <cstdint>
#include <format>
#include <limits>
#include <random>
#include <string>
#include <utility>

#include "gelex/exception.h"

namespace gelex::test
{

bool are_matrices_equal(
    const Eigen::Ref<Eigen::MatrixXd>& mat1,
    const Eigen::Ref<Eigen::MatrixXd>& mat2,
    double tol)
{
    if (mat1.rows() != mat2.rows() || mat1.cols() != mat2.cols())
    {
        return false;
    }

    for (Eigen::Index j = 0; j < mat1.cols(); ++j)
    {
        for (Eigen::Index i = 0; i < mat1.rows(); ++i)
        {
            double val1 = mat1(i, j);
            double val2 = mat2(i, j);

            if (std::isnan(val1) && std::isnan(val2))
            {
                continue;
            }
            if (std::fabs(val1 - val2) > tol)
            {
                return false;
            }
        }
    }
    return true;
}

namespace
{

constexpr std::array<uint8_t, 3> kBedMagicNumber = {0x6C, 0x1B, 0x01};

constexpr std::array<char, 4> kValidNucleotides = {'A', 'C', 'G', 'T'};

constexpr std::array<std::string_view, 24> kChromosomeNames
    = {"1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9",  "10", "11", "12",
       "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X",  "Y"};

uint8_t genotype_to_code(double value)
{
    if (std::isnan(value))
    {
        return 0b01;
    }
    if (value == 0.0)
    {
        return 0b11;
    }
    if (value == 1.0)
    {
        return 0b10;
    }
    if (value == 2.0)
    {
        return 0b00;
    }

    throw ArgumentValidationException(
        std::format(
            "Invalid genotype value: {}, must be 0.0, 1.0, 2.0, or NaN",
            value));
}

}  // namespace

BedFixture::BedFixture() : rng_(std::random_device{}()) {}

std::pair<std::filesystem::path, Eigen::MatrixXd> BedFixture::create_bed_files(
    Eigen::Index num_samples,
    Eigen::Index num_snps,
    double missing_rate,
    double maf_min,
    double maf_max,
    uint64_t seed)
{
    if (num_samples <= 0)
    {
        throw ArgumentValidationException("number of samples must be positive");
    }
    if (num_snps <= 0)
    {
        throw ArgumentValidationException("number of SNPs must be positive");
    }
    if (missing_rate < 0.0 || missing_rate > 1.0)
    {
        throw ArgumentValidationException(
            "missing rate must be in [0.0, 1.0] range");
    }
    if (maf_min < 0.0 || maf_min > 0.5)
    {
        throw ArgumentValidationException(
            "minimum MAF must be in [0.0, 0.5] range");
    }
    if (maf_max < 0.0 || maf_max > 0.5)
    {
        throw ArgumentValidationException(
            "maximum MAF must be in [0.0, 0.5] range");
    }
    if (maf_min > maf_max)
    {
        throw ArgumentValidationException(
            "minimum MAF cannot be greater than maximum MAF");
    }

    rng_.seed(seed);

    std::uniform_real_distribution<double> maf_dist(maf_min, maf_max);
    std::uniform_real_distribution<double> missing_dist(0.0, 1.0);
    std::uniform_real_distribution<double> genotype_dist(0.0, 1.0);

    Eigen::MatrixXd genotypes(num_samples, num_snps);

    for (Eigen::Index snp_idx = 0; snp_idx < num_snps; ++snp_idx)
    {
        double maf = maf_dist(rng_);
        double p = maf * maf;
        double q = 2.0 * maf * (1.0 - maf);
        double r = (1.0 - maf) * (1.0 - maf);

        for (Eigen::Index sample_idx = 0; sample_idx < num_samples;
             ++sample_idx)
        {
            if (missing_dist(rng_) < missing_rate)
            {
                genotypes(sample_idx, snp_idx)
                    = std::numeric_limits<double>::quiet_NaN();
                continue;
            }

            double rand = genotype_dist(rng_);
            if (rand < p)
            {
                genotypes(sample_idx, snp_idx) = 2.0;
            }
            else if (rand < p + q)
            {
                genotypes(sample_idx, snp_idx) = 1.0;
            }
            else
            {
                genotypes(sample_idx, snp_idx) = 0.0;
            }
        }
    }

    auto sample_ids = generate_random_sample_ids(num_samples, rng_);
    auto snp_ids = generate_random_snp_ids(num_snps, rng_);
    auto chromosomes = generate_random_chromosomes(num_snps, rng_);

    current_prefix_ = file_fixture_.generate_random_file_path();
    write_bed_file(genotypes);
    write_bim_file(num_snps, snp_ids, chromosomes);
    write_fam_file(num_samples, sample_ids);

    return {current_prefix_, genotypes};
}

std::pair<std::filesystem::path, Eigen::MatrixXd>
BedFixture::create_deterministic_bed_files(
    const Eigen::MatrixXd& genotypes,
    const std::vector<std::string>& sample_ids,
    const std::vector<std::string>& snp_ids,
    const std::vector<std::string>& chromosomes,
    const std::vector<std::pair<char, char>>& alleles)
{
    const Eigen::Index num_samples = genotypes.rows();
    const Eigen::Index num_snps = genotypes.cols();

    if (!sample_ids.empty()
        && static_cast<Eigen::Index>(sample_ids.size()) != num_samples)
    {
        throw ArgumentValidationException(
            std::format(
                "Sample ID count {} does not match genotype rows {}",
                sample_ids.size(),
                num_samples));
    }
    if (!snp_ids.empty()
        && static_cast<Eigen::Index>(snp_ids.size()) != num_snps)
    {
        throw ArgumentValidationException(
            std::format(
                "SNP ID count {} does not match genotype columns {}",
                snp_ids.size(),
                num_snps));
    }
    if (!chromosomes.empty()
        && static_cast<Eigen::Index>(chromosomes.size()) != num_snps)
    {
        throw ArgumentValidationException(
            std::format(
                "Chromosome count {} does not match genotype columns {}",
                chromosomes.size(),
                num_snps));
    }
    if (!alleles.empty()
        && static_cast<Eigen::Index>(alleles.size()) != num_snps)
    {
        throw ArgumentValidationException(
            std::format(
                "Allele pair count {} does not match genotype columns {}",
                alleles.size(),
                num_snps));
    }

    // 生成默认值（如果未提供）
    std::vector<std::string> final_sample_ids = sample_ids;
    if (final_sample_ids.empty())
    {
        final_sample_ids = generate_random_sample_ids(num_samples, rng_);
    }

    std::vector<std::string> final_snp_ids = snp_ids;
    if (final_snp_ids.empty())
    {
        final_snp_ids = generate_random_snp_ids(num_snps, rng_);
    }

    std::vector<std::string> final_chromosomes = chromosomes;
    if (final_chromosomes.empty())
    {
        final_chromosomes = std::vector<std::string>(num_snps, "1");
    }

    current_prefix_ = file_fixture_.generate_random_file_path();
    write_bed_file(genotypes);
    write_bim_file(num_snps, final_snp_ids, final_chromosomes, alleles);
    write_fam_file(num_samples, final_sample_ids);

    return {current_prefix_, genotypes};
}

std::vector<std::byte> BedFixture::encode_variant(
    const Eigen::VectorXd& variant)
{
    const Eigen::Index n = variant.size();
    const Eigen::Index bytes_per_var = (n + 3) / 4;
    std::vector<std::byte> result(bytes_per_var, std::byte{0});

    for (Eigen::Index i = 0; i < n; ++i)
    {
        uint8_t code = genotype_to_code(variant(i));
        const size_t byte_idx = i / 4;
        const size_t bit_pos = (i % 4) * 2;
        result[byte_idx] |= std::byte(code << bit_pos);
    }

    return result;
}

void BedFixture::write_bed_file(const Eigen::MatrixXd& genotypes)
{
    const Eigen::Index num_samples = genotypes.rows();
    const Eigen::Index num_snps = genotypes.cols();
    const Eigen::Index bytes_per_var = (num_samples + 3) / 4;

    auto bed_path = current_prefix_;
    bed_path.replace_extension(".bed");

    std::vector<std::byte> bed_content;
    bed_content.reserve(3 + (num_snps * bytes_per_var));

    bed_content.push_back(std::byte{kBedMagicNumber[0]});
    bed_content.push_back(std::byte{kBedMagicNumber[1]});
    bed_content.push_back(std::byte{kBedMagicNumber[2]});

    for (Eigen::Index snp_idx = 0; snp_idx < num_snps; ++snp_idx)
    {
        auto encoded = encode_variant(genotypes.col(snp_idx));
        bed_content.insert(bed_content.end(), encoded.begin(), encoded.end());
    }

    [[maybe_unused]] auto bed_file_path
        = file_fixture_.create_named_binary_file(
            bed_path.filename().string(),
            std::span<const std::byte>(bed_content.data(), bed_content.size()));
}

void BedFixture::write_bim_file(
    Eigen::Index num_snps,
    std::span<const std::string> snp_ids,
    std::span<const std::string> chromosomes,
    std::span<const std::pair<char, char>> alleles)
{
    auto bim_path = current_prefix_;
    bim_path.replace_extension(".bim");

    std::string bim_content;
    for (Eigen::Index i = 0; i < num_snps; ++i)
    {
        std::string chrom
            = chromosomes.empty() ? "1" : std::string(chromosomes[i]);
        std::string snp_id
            = snp_ids.empty() ? std::format("rs{}", i + 1) : snp_ids[i];

        std::pair<char, char> allele_pair;
        if (alleles.empty())
        {
            allele_pair = generate_random_alleles(rng_);
        }
        else
        {
            allele_pair = alleles[i];
        }

        bim_content += std::format(
            "{} {} {} {} {} {}\n",
            chrom,
            snp_id,
            static_cast<double>((i + 1) * 100),
            i + 1,
            allele_pair.first,
            allele_pair.second);
    }

    [[maybe_unused]] auto bim_file_path = file_fixture_.create_named_text_file(
        bim_path.filename().string(), bim_content);
}

void BedFixture::write_fam_file(
    Eigen::Index num_samples,
    std::span<const std::string> sample_ids)
{
    auto fam_path = current_prefix_;
    fam_path.replace_extension(".fam");

    std::string fam_content;
    for (Eigen::Index i = 0; i < num_samples; ++i)
    {
        std::string sample_id = sample_ids.empty()
                                    ? std::format("sample{}", i + 1)
                                    : sample_ids[i];

        std::string fid = std::format("fam{}", (i % 5) + 1);

        fam_content
            += std::format("{} {} 0 0 {} -9\n", fid, sample_id, (i % 2) + 1);
    }

    [[maybe_unused]] auto fam_file_path = file_fixture_.create_named_text_file(
        fam_path.filename().string(), fam_content);
}

std::pair<char, char> BedFixture::generate_random_alleles(std::mt19937_64& rng)
{
    std::uniform_int_distribution<size_t> dist(0, kValidNucleotides.size() - 1);

    char a1 = kValidNucleotides.at(dist(rng));
    char a2 = kValidNucleotides.at(dist(rng));

    while (a1 == a2)
    {
        a2 = kValidNucleotides.at(dist(rng));
    }

    return {a1, a2};
}

std::vector<std::string> BedFixture::generate_random_sample_ids(
    Eigen::Index num_samples,
    std::mt19937_64& rng)
{
    std::vector<std::string> ids;
    ids.reserve(num_samples);

    for (Eigen::Index i = 0; i < num_samples; ++i)
    {
        ids.push_back(std::format("sample{}", i + 1));
    }

    return ids;
}

std::vector<std::string> BedFixture::generate_random_snp_ids(
    Eigen::Index num_snps,
    std::mt19937_64& rng)
{
    std::vector<std::string> ids;
    ids.reserve(num_snps);

    for (Eigen::Index i = 0; i < num_snps; ++i)
    {
        ids.push_back(std::format("rs{}", i + 1));
    }

    return ids;
}

std::vector<std::string> BedFixture::generate_random_chromosomes(
    Eigen::Index num_snps,
    std::mt19937_64& rng)
{
    std::vector<std::string> chromosomes;
    chromosomes.reserve(num_snps);

    std::uniform_int_distribution<size_t> dist(0, kChromosomeNames.size() - 1);

    for (Eigen::Index i = 0; i < num_snps; ++i)
    {
        chromosomes.emplace_back(kChromosomeNames.at(dist(rng)));
    }

    return chromosomes;
}

}  // namespace gelex::test
