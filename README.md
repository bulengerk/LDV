
# LDV

Ligo Data Visualizer


## Where to find gravitonal events?

[List of observed events](https://gwosc.org/eventapi/html/allevents/)


## Usage

### Build
```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake libhdf5-dev libcurl4-openssl-dev
mkdir build && cd build
cmake ..
cmake --build . -j
```


### Run app
Choose your gravitational data from gwosc and call:

```bash
  ./gwosc_plotter https://gwosc.org/archive/data/DiscO4a_16KHZ_R1/1384120320/H-H1_GWOSC_DiscO4a_16KHZ_R1-1384779776-4096.hdf5
```




e.g.

check the output of the out directory in order to find your visualization


