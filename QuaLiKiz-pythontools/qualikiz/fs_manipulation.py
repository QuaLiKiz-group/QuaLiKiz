"""
Copyright Dutch Institute for Fundamental Energy Research (2016)
Contributors: Karel van de Plassche (karelvandeplassche@gmail.com)
License: CeCILL v2.1
"""
import os

from .qualikizrun import QuaLiKizBatch


def move_completed(searchdir, targetdir):
    batchlist = QuaLiKizBatch.from_dir_recursive(searchdir)
    for batch in batchlist:
        if batch.is_done():
            batchdir = os.path.join(batch.batchsdir, batch.name)
            os.rename(batchdir, os.path.join(targetdir, batch.name))
