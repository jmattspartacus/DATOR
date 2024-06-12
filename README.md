\mainpage DATOR: Data Analysis Toolkit for ORRUBA

This is a toolkit of software designed to assist with data analysis for ORRUBA experiments. Specifically, it aims to handle cases where ORRUBA is coupled to additional detector systems such as Gretina, the S800, or others.

In practice, what this means is the toolkit can handle data from *any* combination of sources (not necessarily including ORRUBA), as long as that data is packaged in the well-established "Global Event Builder", or GEB format. However, since the toolkit has been developed by the ORRUBA collaboration with a view to ORRUBA analysis as a first use-case, the name remains.

# Introduction to DATOR

## Data pipeline

The data pipeline for DATOR is shown below.

\image html DATOR_pipeline.png
\image latex DATOR_pipeline.png

### Creating the merged Global.dat

The raw data comes in the form of several files, typically a Global.dat or HFC.dat from one or more auxiliary systems, and a *.ldf file from the ORRUBA/Orphas data acquisition. A standalone executable, LDFMerge is included in DATOR, which merges the Global.dat with a single *.ldf file. This assumes that both files are already well time-ordered, and goes through event-by-event, removing any Type-19 (ORRUBA) events from the Global.dat, and replacing them with new Type-19 events constructed from the *.ldf data.

If additional files or formats are used in the future (i.e. some auxiliary system which uses neither the Orphas LDF nor is merged into the Global.dat), additional programs will be needed to complete the merge of these data streams.

None-the-less, the end product of this stage is a "GlobalMerged.dat", which contains all the data in a single GEB file, time-ordered, and ready for physics analysis. This merging requires no decisions to be made, and should be done only once. After the GlobalMerged.dat is created, this can be treated as the definitive version of the data and all analysis can proceed on this file.

### Time correlation

The next step of the analysis is denoted "time correlation", and may also be referred to as "event building". Here, the GlobalMerged.dat is traversed with a rolling coincidence window, to group sub-events which are close together in time into one "Physics event". This takes place in the ```Reader::Next()``` function, is a fairly standard procedure, and shouldn't need to be altered much if at all in a typical analysis. The coincidence window can be specified by changing the ```Reader::coinc_window``` variable.

### Physics processing

The real magic of analysis happens in the "physics processing" stage. Here each GEB sub-event is read in, and passed to a user-defined "processor" object, according to its GEB type. So for example, a Type-1 event (Gretina signal decomp) would be passed to a Type-1 processor that estimates the initial interaction point for the Gretina hit, and transforms that into a global coordinate system. A Type-19 event (ORRUBA data) could be processed by conducting a pedestal subtraction and gain alignment, matching front- and back- hits from the same detector, and calculating the position of the detected particle in the global coordinate system. All the derived information and associated methods for each type of processing are stored within the processor object and can be subsequently accessed for further processing and/or outputs.

\image html DATOR_physics.png
\image latex DATOR_physics.png

The processor objects must have three methods defined:

``` Processor::Process(unsigned long long int timestamp, unsigned short int *data, unsigned short int length); ```

This takes the payload of the GEB event (pointed to by *data), and does the processsing.

```Processor::ProcessFinal();```

This is any additional processing to take place at the end of the physics event, after all subevents have been read and processed. For example, add-back for Gretina.

```Processor::Reset();```

This is called at the beginning of each physics event: i.e. delete or reset any data from the previous event. Optionally, a fourth member can be defined

```Processor::PrintSummary(std::ostream &out);```

which will be called at the end of each file processed. This is useful to print a summary of statistics and diagnostics to ```std::cout``` or a log file.

### The sort program

The time correlation, physics processing, and eventual output to either histograms or other event-by-event data formats (i.e. ROOT Trees) all takes place in a user-compile executable, i.e. the "sort code". A block diagram showing a typical case is shown below.

\image html DATOR_sort.png
\image latex DATOR_sort.png

The best way to understand this in detail is to work through an example. Take a look at BasicSort.cc and follow the instructions there to compile and run the sort program on your own system. You can use it as a template to add your own output histograms, change the diagnostic reporting, and swap in and out processors for your own specific needs.

### Writing your own processor

Several standard procesors for common data types are included in DATOR, however it is likely that at some point you will need to alter or even write your own from scratch. A tutorial to get you started on this can be found here (see MyProcessor.cc).

### The pruned Global.dat

A useful feature provided by the DATOR::Reader class is the ability to write out an entire physics event to a new file, in the original GEB/Global.dat format. The new file is referred to as the "pruned" Global.dat. A typical use case for this feature would be to set some cuts or conditions on what qualifies a "good" physics event. Perhaps this relates to an incoming beam PID gate, or a multiplicity condition above some energy threshold, or some coincidence between several different detector systems. In any case, it may be the case that only a small fraction of the events in the raw data file actually contain useful or relevant information. Thus, after basic processing is complete, we can reject all other events and only write out those that are of interest to a separate file. 

Since it is in the original format, it can be processed through exactly the same sort program, to give exactly the same output histograms. However, if its size is significantly smaller, it will be significantly faster to process as a result. Thus, one can set rough cuts and conditions to isolate the interesting events, and then fine-tune the analysis on the subset of the data quickly. There are several parameters of the DATOR::Reader class that relate to the file name and location of the pruned Global.dat, and then writing the physics events to this file is achieved by calling DATOR::Reader::Write() while inside the event loop.

# Installation

Clone the repository into a new folder. Once inside the repository, type

~~~~~~~~~~~~~~{.sh}
make
make install
~~~~~~~~~~~~~~

To build and install DATOR. There are a few details inside the Makefile which may be of note:

    INSTALLDIR=$(HOME)/.local

This line sets where the built libraries, binaries, and header files will be copied to. The libraries go into

    $(INSTALLDIR)/lib/

The binaries into

    $(INSTALLDIR)/bin/

and the headers into

    $(INSTALLDIR)/include/DATOR/

I use ```$(HOME)/.local/``` as a personal preference, but you can choose anywhere on your system. Just make sure that in your ```.bashrc``` file (or equivalent), you have the binary folder in your ```$PATH``` variable, and the library folder in your ```$LD_LIBRARY_PATH``` (Linux) or ```$DYLD_LIBRARY_PATH``` (macOS).

The ```CFLAGS``` variable defined in the Makefile contains a reference to ```root-config --cflags``` which will ensure DATOR is built with the same flags as your local ROOT installation. However, if you don't have or need ROOT, just remove this part from the definition and you should be fine. 

# Contributors

- Tim Gray (University of Tennesse, Knoxville/ORNL, tgray30@utk.edu)
