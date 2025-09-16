import time
from functools import wraps

from gelexy._logger import logger


def timeit(N: int = 10) -> None:
    """
    Decorator to measure and log the execution time of a function over multiple runs.

    :param N: Number of times to execute the decorated function. Defaults to 10.
    :return: The decorator that wraps the function, logs timing info, and returns the result of the last call.

    Usage:
        @timeit(N=5)
        def my_function():
            # function body
            pass

        my_function()
    """

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            logger.info("--- Running '%s' %d times ---", func.__name__, N)

            start_time = time.perf_counter()

            result = None
            for _ in range(N):
                result = func(*args, **kwargs)
            end_time = time.perf_counter()

            total_time = end_time - start_time
            avg_time = total_time / N

            logger.info("Total time %.6f s", total_time)
            logger.info("Average time %.6f s", avg_time)
            logger.info("--- Finished ---")
            return result

        return wrapper

    return decorator
