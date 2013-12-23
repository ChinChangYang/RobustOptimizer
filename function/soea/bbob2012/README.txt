This is the Matlab version of the COCO software
    --- http://coco.gforge.inria.fr/ ---

It was tested using with Matlab r2008b and Octave

As usual, this software comes with absolutely no warranty (see Warranty at end).

It is here assumed that you know what COCO/BBOB is about http://coco.gforge.inria.fr/doku.php?id=bbob-2010 and that you have read the corresponding documentation http://coco.gforge.inria.fr/

Installation
============
If you are reading this, you already succeeded in unpacking the COCO tar file :-)

After untarring the archive, the bbob.vXX.XX/matlab directory should contain the following files:

benchmarkinfos.txt
benchmarks.m
benchmarksnoisy.m
exampleexperiment.m
exampletiming.m
fgeneric.m
MY_OPTIMIZER.m
README.txt

* exampleexperiment.m - that launches a very quick but complete experiment using the random MY_OPTIMIZER optimizer (see code).

* exampletiming.m - that launches the timing experiment (warning, this program takes at least half a minute before the first output to console, be patient :-)

If some run-time errors appear when running either of those programs, get in touch immediately with bbob@lri.fr indicating all information about your system, as well as the complete error message(s).

Customization
=============
you might at some point want to try your own optimizer: this is fairly easy.

1- edit exampleexperiment.m and replace the strings "PUT_..." with the correct names. If you modify the 'PUT_MY_BBOB_DATA_PATH' string, you need to run the createfolder.py python script with the corresponding new string, see above.

2- write a MY_OPTIMIZER.m function that sets up everything for your optimizer, and calls it with the correct arguments

You're done :-)

In case of problem, send a mail to bbob@lri.fr but please make sure first that your optimizer runs smoothly as a standalone application.
To discuss more general issues about BBOB: bbob-discuss@lists.lri.fr

Acknowledgments
===============
The BBOBies would like to acknowledge Miguel Nicolau who contributed by testing this version of the code.
We would also like to acknowledge Alvaro Fialho who helped transposing the C-version to the C++ language.

Warranty
========
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL ANY CONTRIBUTOR TO THIS SOFTWARE BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

