#include "chenx/predictor.h"

#include <string_view>

#include <armadillo>

/*namespace chenx*/
/*{*/
/*using arma::rowvec;*/
/**/
/*Predictor::Predictor(*/
/*    std::string_view train_bed_file,*/
/*    std::string_view method,*/
/*    rowvec&& center,*/
/*    double scale_factor,*/
/*    uint64_t chunk_size)*/
/*{*/
/*    if (method == "add")*/
/*    {*/
/*        cross_grm_ = std::make_unique<AddCrossChunkGrm>(*/
/*            train_bed_file, std::move(center), scale_factor, chunk_size);*/
/*    }*/
/*    else if (method == "dom")*/
/*    {*/
/*        cross_grm_ = std::make_unique<DomCrossChunkGrm>(*/
/*            train_bed_file, std::move(center), scale_factor, chunk_size);*/
/*    }*/
/*    else*/
/*    {*/
/*        throw std::invalid_argument("Invalid method");*/
/*    }*/
/*}*/
/**/
/*Predictor::Predictor(*/
/*    std::string_view train_bed_file,*/
/*    std::string_view method,*/
/*    rowvec&& center,*/
/*    double scale_factor)*/
/*{*/
/*    if (method == "add")*/
/*    {*/
/*        cross_grm_ = std::make_unique<AddCrossGrm>(*/
/*            train_bed_file, std::move(center), scale_factor);*/
/*    }*/
/*    else if (method == "dom")*/
/*    {*/
/*        cross_grm_ = std::make_unique<DomCrossGrm>(*/
/*            train_bed_file, std::move(center), scale_factor);*/
/*    }*/
/*    else*/
/*    {*/
/*        throw std::invalid_argument("Invalid method");*/
/*    }*/
/*}*/
/**/
/*}  // namespace chenx*/
