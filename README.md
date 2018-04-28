# b26_toolkit
b26_toolkit is an implementation of the scripts
and instruments used by the NV Nanomechanics project in the Lukin Group 
at Harvard University. It contains PyLabControl-compatible instruments 
and scripts that are utilized in the project, as well well as plotting
and data processing convenience functions.

b26_toolkit is built on python 3.6.x. It was built by Arthur Safira,
Jan Gieseler, and Aaron Kabcenell in the Lukin Group at Harvard University. 
It is distributed under the GPLv3 license. For more information, see LICENSE.txt .


## Getting Started

### Installation
The simplest way to install b26_toolkit is with the command-line utility pip. To install simply issue the command

```>>> pip install b26_toolkit```

in the windows command-line interface.

#### Exporting Instruments and Scripts for use in PyLabControl
In order to use the implemented instruments and scripts, one must first export them into appropriate .b26 files . To export them, use

``` >>> PyLabControl --export <file_path> ```

in the windows commandline interface, where `<file_path>` is a valid path to save the exported scripts and instruments.
The instruments and scripts that were able to be exported will then be located at that file path, and could be imported into the PyLabControl GUI.


## FAQ
+ **Some instruments seemed to have exported; How can I see what errors occured?**

To see the errors that occur when an instrument does not export, simply try importing that instrument and creating an instance of it.
For example, while in a python interpreter, one could write

```
import b26_toolkit.PiezoController
piezo_controller_x = PiezoController('x')
```
and see what errors occur. Feel free to reach out to the authors with any difficulties.

+ **What instruments have been implemented?**

For a full list of instruments, see the files in b26_toolkit/src/instruments/ .

+ **How can we send feedback?**

Feel free to create an issue on the issue tracker if you find any bugs, or contact the authors directly.
