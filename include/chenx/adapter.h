#include <armadillo>
namespace chenx
{
template <typename From, typename To>
struct Adapter
{
    static To to(From from)
    {
        static_assert(sizeof(From) == 0, "Adapter not implemented");
    }
};

template <typename T>
struct Adapter<T, T>
{
    static T&& to(T&& value) noexcept { return std::forward<T>(value); }
};
}  // namespace chenx
