"""Utility functions for the gelexy package."""

from pathlib import Path


def valid_path(
    input_path: str | Path, suffixes: tuple[str] | None = None
) -> Path:
    """
    Validate the input path and checks for specified suffixes.

    Parameters
    ----------
    input_path : str | Path
        The path to validate. Can be a string or a Path object.
    suffixes : tuple[str] | None, optional
        A tuple of valid suffixes. If provided, the function checks if the path has one of these suffixes.

    Returns
    -------
    Path : str
        the valid full path.

    Raises
    ------
    TypeError
        If the input is not a string or a Path object.
    FileNotFoundError
        If the path does not exist.
    ValueError
        If the path is neither a file nor a directory, or if the suffix is not valid.
    """
    path = Path(input_path) if isinstance(input_path, str | Path) else None

    if path is None:
        msg = "Input must be a string or a Path object."
        raise TypeError(msg)

    if not path.exists():
        msg = f"The path '{path}' does not exist."
        raise FileNotFoundError(msg)

    if not (path.is_file() or path.is_dir()):
        msg = f"The path '{path}' is neither a file nor a directory."
        raise ValueError(msg)

    if suffixes:
        suffixes = [s if s.startswith(".") else f".{s}" for s in suffixes]
        if path.suffix not in suffixes:
            msg = f"The file '{path}' does not have a valid suffix. Expected one of: {suffixes}"
            raise ValueError(msg)
    return str(path.absolute())
