# Hetero 

### (Under development)
                     
              
### Install:
When installing from GitHub source the following tools are required:

* [Autoconf](http://www.gnu.org/software/autoconf)
* [Automake](http://www.gnu.org/software/automake)

To generate the configure script and make files:

```
./autogen.sh
```

To compile and install hetero in /usr/local:

```
$ ./configure
$ make
$ sudo make install
```
To install hetero in a specified directory:

```
$ ./configure --prefix=/opt/hetero
$ make 
$ make install 
```

Hetero uses OpenMP for parallelization, which requires a modern compiler such as GCC 4.2 or greater. If you have an older compiler, it is best to upgrade your compiler if possible. If you have multiple versions of GCC installed, you can specify a different compiler:

```
$ ./configure CC=gcc-xx CXX=g++-xx 
```

For the best performance of hetero, pass `-O3` flag:  

```
$ ./configure CFLAGS='-g -O3' CXXFLAGS='-g -O3' 
```

To run hetero, its executables, `CreateBloom` and `hetero`, should be found in your PATH. If you installed hetero in /opt/hetero, add /opt/hetero/bin to your PATH:

```
$ PATH=/opt/hetero/bin:$PATH
```

