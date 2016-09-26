import os

from .qualikizrun import QuaLiKizRun, QuaLiKizBatch

def move_completed(searchdir, targetdir):
    batchlist = QuaLiKizBatch.from_dir_recursive(searchdir)
    for batch in batchlist:
        if batch.is_done():
            batchdir = os.path.join(batch.batchsdir, batch.name)
            os.rename(batchdir, os.path.join(targetdir, batch.name))
