# climate-companions-bipartite
Code for investigating the event co-participation based bipartite graph structure in the Ilmastokumppanit (Climate Companions) city-to-business network. Should be applicable also to other bipartite module analyses.

<code>functions.py</code> contains functions needed for creating and analysing the bipartite network. The functions can be called from either a frontend script (<code>fronend.py</code> was used for creating the results of the related article) or interactively. <code>parameters.py</code> stores all parameters used in the anlysis.

For a short(-ish) walkthrough of the code and main results, have a look at <code>walkthrough.py</code> (related data is stored at <code>walkthrough_data</code>). The easiest way to do that is to launch Binder (that saves you the trouble of setting package versions or even installing python):

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/onerva-korhonen/climate-companions-bipartite/HEAD)

(Preparing the environment may take a while. When it opens, just click <code>walkthrough.py</code> to launch the notebook).

The versions used are Python 2.7, Pandas 0.22.0, Numpy 1.14.0, NetworkX 2.1, Scipy 1.0.0, and Matplotlib 2.1.2.

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.
