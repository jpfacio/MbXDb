import shutil
from pathlib import Path


def clear_working(
    outdir: Path = Path("tmp/ncbi_tmps"),
    subdirs: tuple = ("srr", "assembly", "bins"),
) -> None:
    """Remove the per-chunk working directories, keeping the directory layout.

    Used between chunks to reclaim disk after each batch of SRRs has been
    processed and its bins moved out. Every listed subdirectory is deleted and
    recreated empty so the next chunk starts from a clean slate.

    Args:
        outdir (Path): root of the NCBI working directories
        subdirs (tuple): subdirectory names to clear

    Returns:
        None
    """

    outdir = Path(outdir)

    for name in subdirs:
        target = outdir / name
        if target.exists():
            shutil.rmtree(target)
        target.mkdir(parents=True, exist_ok=True)
