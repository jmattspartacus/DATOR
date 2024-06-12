/*!
  \file MyProcessor.cc

  This example file will help you to make your own processor for a GEB data type. As an example, we will use Type-19 ORRUBA data. We will start by making two files, a header file "MyProcessor.hh" and an implementation file "MyProcessor.cc". Put these in the same location as your working source directory, wherever BasicSort.cc or your main sort program source is. Inside the header file, we define a new class which inherits from the DATOR::Processor class, defined in Reader.hh. 

~~~~~~~~~~~~~~~{.cpp}
class MyProcessor : public DATOR::Processor {
  public:
  void Process(unsigned long long int timestamp, unsigned short int *data, unsigned short int length);
  void ProcessFinal();
  void Reset();
  void PrintSummary(std::ostream &out);

  unsigned long long int time;
  std::vector<double> chans;
  std::vector<double> vals; 

  //here go the derived data things
  int threshold;
  void MakeSum();
  double totalSum;
};
~~~~~~~~~~~~~~~

I've also defined a member variable ```totalSum```, and a member function ```MakeSum()``` which are stand-ins for any kind of derived data structure and processing function. The members could be other classes, collections (i.e. vector), or really anything at all. You can have as many of these as you like to represent whatever complexity of data you need. Next, let's define the implentation in MyProcessor.cc:

~~~~~~~~~~~~~~~{.cpp}
#include "MyProcessor.hh"

void MyProcessor::Process(unsigned long long int timestamp, unsigned short int *data, unsigned short int length) {
  //here goes the code to take the data from the C-style vector *data (length bytes), and import it into the vectors chans and vals
  //I'll leave you to fill this in, if you get stuck look through Reader.cc and the Orruba::Basic class.
  //The type-19 data is structured as a sequence of channel-value pairs, each as an unsigned short int. The last two are a "footer", and should always be the same: (0xFFFF,0xFFFF).
}

void MyProcessor::ProcessFinal() {
  //as a silly example, I'll say that the processing step we wish to take is sum up all the vals for channels greater than some threshold, and store these in totalSum. 
  //let's do this in the MakeSum() function though, and simply call that function from here.

  MakeSum();
}

void MyProcessor::Reset() {
  //when we reset the event we want to clear the two vectors, and reset totalSum to be zero
  chans.clear();
  vals.clear();
  totalSum = 0;
}

void MyProcessor::MakeSum() {
  //let's do our summing here. Remember this is a pretty trivial example, but this function could have any amount of complexity, calling other member functions, and so on an so forth.

  for (int i=0; i<chans.size(); ++i) {
    if (chans[i] >= threshold) { 
      totalSum += vals[i];
    }
  }
}
~~~~~~~~~~~~~~~

That wraps up the implementation of our processor class. Now we need to actually put the processor in the sort code. So switching to BasicSort.cc (or your sort source file), we need to create the processor:

~~~~~~~~~~~~~~~{.cpp}
MyProcessor myproc;
myproc.threshold = 50; //this is the threshold for the sum. In a more realistic example, this might be any sort of configuration variable: a channel map, geometry definitions, gain calibrations
~~~~~~~~~~~~~~~

Then we can add it to the reader, associated with type-19 events
~~~~~~~~~~~~~~{.cpp}
reader.AddProcessor(19, &myproc);
~~~~~~~~~~~~~~

And now the ```MyProcessor::Process()``` and ``` MyProcessor::ProcessFinal()``` functions we defined above will be called automatically during the ```reader.Next()``` event loop. Now we can do something with the results. So inside the event loop, add a line

~~~~~~~~~~~~~~{.cpp}
std::cout << myproc.totalSum << std::endl;
~~~~~~~~~~~~~~

Which will spam the screen with numbers, but if they're all different and non-zero, that indicates we're summing correctly for each event. For bonus points, you can add a separate member variable which counts the number of channels summed and print that out as well, or make a histogram of the summed values.

To get this to compile we'll need to make a small change to the Makefile as well:

    BasicSort : BasicSort.cc libReader libGRETINA libORRUBA
    	$(CC) $(CFLAGS) -I$(INSTALLDIR)/include/DATOR/ -o BasicDatorSort BasicSort.cc MyProcessor.cc $(ROOTLIBS) $(LIBS) -L$(INSTALLDIR)/lib -lReader -lGRETINA -lORRUBA

simply adding the implementation file MyProcessor.cc to the source files passed to the compiler. You should be able to make and run BasicSort as usual.

*/

void MyProcessor::Process(unsigned long long int timestamp, unsigned short int *data, unsigned short int length) {};
void MyProcessor::ProcessFinal() {};
