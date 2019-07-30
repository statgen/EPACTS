# EPACTS - Efficient and Parallelizable Association Container Toolbox

### What is EPACTS?

EPACTS is a versatile software pipeline to perform various statistical tests for identifying genome-wide association from sequence data through a user-friendly interface, both to scientific analysts and to method developer.s

### Downloading and Installing EPACTS

Prerequisite Packages
- zlib
- ghostscript
- R-studio
- groff
- gnuplot
- automake  

Installing these packages on Centos:
<pre>
$ sudo yum install epel-release
$ sudo yum update
$ sudo yum install zlib ghoshcript R groff gnuplot automake -y
</pre>

Initialize the autoreconf program:
<pre>
$ autoreconf -f -i
</pre>

You can clone the current snapshot of this repository to install as well

<pre>
$ git clone https://github.com/statgen/EPACTS.git
$ cd EPACTS
$ ./configure --prefix [/path/to/install]
$ make
$ make install
</pre>

### EPACTS Documentation

The latest version of documentation of EPACTS can be found at
http://genome.sph.umich.edu/wiki/EPACTS

### Feedbacks

Feel free to contact Hyun Min Kang (hmkang@umich.edu) or joint EPACTS Google group (http://groups.google.com/group/epacts) to ask questions about EPACTS
