#ifndef COMMANDLINEPARSER_H_
#define COMMANDLINEPARSER_H_

#include <string>

class CommandLineParser {
public:
    static std::string parse(int argc, char **argv);
};

#endif
