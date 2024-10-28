# LT-SPRAT quick reduction pipeline

Performs bulk quick extraction of LT-SPRAT spectra using the 1D output from the LT pipeline. A flux correction is applied depending on airmass.

Alternatively, the notebook provides a way for more manual extracting if desired.

Requires Pyraf, see Astroconda

## Initial setup

To run, the code just needs to be cloned or downloaded.

Download the required .fits files from the LT data archive and save in the folder *InputSpectra*.

If the automatic LT reduction has found a good trace, QuickSPRATPipeline.py will quickly extract it.

If you want to extract another trace (e.g. to inspect the host that happened to be in the same slit), use the first part of the notebook. This can of course also be used instead of the QuickSPRATPipeline.py

If the LT reduction found no trace (the file does not end with 2.fits, but 1.fits) the second part of the notebook can be used for trace extraction. Note that this also requires the arc (usually the next file, the second letter should be an 'a' instead of 'e').

## Commands

In an IRAF-ready conda environment, navigate to the LT-Quick-Bulk folder and type 

```python QuickSPRATPipeline.py```

There are also other tags that can be applied, but remember that as this script acts on all the files these functions will be applied to every spectrum reduced.

To specify a redshift to plot at, use
```python QuickSPRATPipeline.py -z <redshift>```

To apply an extinction correction during plotting:
```python QuickSPRATPipeline.py -E <E(B-V)>```

Now it is possible to do a basic cosmic ray removal on the final file. Be careful with this however because we are dealing with 1D spectra it will also remove host galaxy lines as it cannot tell the difference

```python QuickSPRATPipeline.py -C```

## Output

The output is an initial extraction in *InputSpectra* and two plots in *plots* which compare the initial and final reductions.

The final output is a text file with the extension **\*.w.txt** in the folder _OutputSpectra_