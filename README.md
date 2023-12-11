# GenUtahSky
Radiance implementation of Preetham, Shirley, Smits model for sky color, plus more

### What is this?

Users of the [Radiance Synthetic Imaging System](http://radsite.lbl.gov/radiance/) can use the software in this repository to generate and render realistic skies for anywhere on Earth at any requested time (no clouds, sorry). The colors of the sky come from the Preetham, Shirley, Smits (1999) paper. Extras added by this author include: realistic sun disk size and color, moon with proper brightness and position (but not phase), Jupiter and Venus, plus a glowing starfield at night.

Here are [slides from a presentation about Radiance skies](http://markjstock.org/radiance/radiance_harvard_09_stock.pdf) and a [video of 24 hours of skies over the Grand Canyon](https://www.youtube.com/watch?v=BJWviVcu_Qo).

### How to use it

To build the genutahsky executable, you will need a C compiler and the libnova libraries. On Fedora or Ubuntu, you can get this by running, respectively:

    sudo dnf install libnova-devel
    sudo apt-get install libnova-dev

Once you have those, grab this repository with:

    git clone --depth 1 https://github.com/markstock/GenUtahSky.git
    cd GenUtahSky
    make

Run a simple test with

    genutahsky 6 23 19

To use the system in a Radiance scene, simply copy all of the files in this directory into the directory containing your scene, and add a command like the following to one of your `.rad` files:

    !genutahsky 6 23 19 -t 5.0 -a 45 -o 105

This will generate the proper materials for a sky at 45N 105W, on June 23 at 7pm local time, with a turbidity of 5 (2 represents very clear air, 10 is hazy).

If you don't want to clutter your working directory up, drop these files into your Radiance `lib` directory, often at `/usr/local/lib/ray`.

If you're going to use `pcond` on night scenes, use `-v -s` instead of `-h`. The latter will blur the stars unrealistically.

### To Do

* Use CMake to build and install the package (putting necessary files into radiance library location and binaries into /usr/local/bin
* Update starfield to https://svs.gsfc.nasa.gov/4851/
* Have build script retrieve the starfield and reduce its resolution (and increase bit depth) instead of storing it directly in the repository

### Copyrights and licensing

The genutahsky.c program uses code from [Radiance](http://radsite.lbl.gov/radiance/), (c) Greg Ward Larson.

All sky color numbers are from the original [Preetham, Shirley, Smits paper](http://www.cs.utah.edu/~shirley/papers/sunsky/sunsky.pdf).

The program optionally links with the [LibNova library](http://libnova.sourceforge.net/) for astronomical calculations, maintained by Liam Girdwood and Petr Kubanek. It also uses [Astronomy](https://github.com/cosinekitty/astronomy) if libnova is unavailable.

The starfield is courtesy the NASA/Goddard Space Flight Center Scientific Visualization Studio, and was generated from the Tycho II Star Survey. I processed their original high-resolution 8bpp image to a smaller 12bpp hdr image. [Original data](http://svs.gsfc.nasa.gov/cgi-bin/details.cgi?aid=3572)

All other material is (c) 2009,23 Mark J. Stock

