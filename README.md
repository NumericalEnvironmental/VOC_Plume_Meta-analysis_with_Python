# VOC_Plume_Meta-analysis_with_Python

![Preview](https://numericalenvironmental.files.wordpress.com/2017/04/cal_general_map.jpeg?w=655&h=463)

This is a python script I have written to read and parse a large groundwater quality database (thousands of wells, with hundreds of thousands of sample events) into individual groundwater plumes, subject to subsequent spatial analyses to develop a simple set of plume metrics. A more detailed description and some example results are provided in my blog:

https://numericalenvironmental.wordpress.com/2017/04/17/meta-analysis-of-over-1000-groundwater-chlorinated-hydrocarbon-plumes-using-python-tools/

The required python 2.7-compatible libraries required to run this script include:
* numpy
* pandas
* scikit-learn
* pyproj
* scipy
* seaborn

In addition, the following input files are needed:

* gama_df.csv: a comma-delimited text file which is an already-produced extraction of chlorinated volatile organic compound (CVOC) data from California’s Groundwater Ambient Monitoring and Assessment (GAMA) database, if this option is selected when the code is run. Otherwise, tab-delimited text files for complete data from all 50+ California counties is required (several GB worth of data, so not provided here). If this latter option is selected, the gama_df.csv file will be created automatically for subsequent runs.

* analytes.txt: a text file listing the CVOCs of interest to read from the database (naming must be consistent with what is used in GAMA).

* counties.txt: a text file listing the names of individual tab-delimited text files corresponding to California counties, as represented in GAMA. This file is required regardless of whether data are read from gama_df.csv summary file or not.

Email me with questions at walt.mcnab@gmail.com.

THIS CODE/SOFTWARE IS PROVIDED IN SOURCE OR BINARY FORM "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

