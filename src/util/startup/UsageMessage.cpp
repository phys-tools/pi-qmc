#include "UsageMessage.h"
#include <iostream>

void UsageMessage::print() {
    std::cout << "usage: pi [OPTIONS] [pimc.xml]" << std::endl;
    std::cout << "  OPTIONS" << std::endl;
    std::cout << "     -h, --help           Print a usage message and exit"
            << std::endl;
    std::cout << "     -V, --version        Print version number and exit"
            << std::endl;
    std::cout << "     -d, --demo=name      Write input file for a demo\n";
    std::cout << "                          (omit name to get a list of "
            << "available demos)" << std::endl;
}
