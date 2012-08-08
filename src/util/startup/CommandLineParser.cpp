#include "CommandLineParser.h"

#include "base/Help.h"
#include "demo/Demo.h"
#include "util/startup/UsageMessage.h"

#if HAVE_GETOPT_H
#include <getopt.h>
#else
#include "util/gnugetopt.h"
#endif

#include <cstdlib>

std::string CommandLineParser::parse(int argc, char** argv) {
    // Parse command line.
    static struct option longopts[] = { { "help", no_argument, NULL, 'h' }, {
            "demo", optional_argument, NULL, 'd' }, { "version", no_argument,
            NULL, 'V' }, { NULL, 0, NULL, 0 } };
    char ch;
    while ((ch = getopt_long(argc, argv, "hdV", longopts, NULL)) != -1) {
        std::string demoName;
        switch (ch) {
        case 'd': {
            if (optarg != NULL) {
                demoName = std::string(optarg);
                std::cout << "Requested demo: " << demoName << "\n"
                        << std::endl;
                Demo *demo = Demo::getDemo(demoName);
                if (demo) {
                    demo->generate();
                    exit(0);
                } else {
                    Demo::listDemos(std::cout);
                    exit(-1);
                }
            } else {
                Demo::listDemos(std::cout);
                exit(0);
            }
        }
            break;
        case 'V':
#ifdef CONFIG_FLAGS
            std::cout<<"Compiled with options: "<<std::endl<<CONFIG_FLAGS<<std::endl<<std::endl;
#endif
            exit(-1);
            break;
        case 'h':
            Help::printHelp();
            break;
        default:
            UsageMessage::print();
            exit(0);
            break;
        }
    }
    argc -= optind;
    argv += optind;
    std::string xmlFileName = "pimc.xml";
    if (argc == 1)
        xmlFileName = argv[0];

    return xmlFileName;
}
