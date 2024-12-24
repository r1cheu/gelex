#include <cmath>
#include <iostream>

class RemlSummary
{
  private:
    double logLikelihood; // 对数似然值
    int numParams;        // 参数个数
    int numSamples;       // 样本量

  public:
    // 构造函数
    RemlSummary(double logL, int params, int samples)
        : logLikelihood(logL), numParams(params), numSamples(samples)
    {
    }

    // 计算AIC
    double computeAIC() const
    {
        return 2 * numParams - 2 * logLikelihood;
    }

    // 计算BIC
    double computeBIC() const
    {
        return std::log(numSamples) * numParams - 2 * logLikelihood;
    }

    // 输出AIC和BIC
    void printSummary() const
    {
        std::cout << "AIC: " << computeAIC() << std::endl;
        std::cout << "BIC: " << computeBIC() << std::endl;
    }
};
