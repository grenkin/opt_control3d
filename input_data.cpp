#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <conio.h>
#include "input_data.h"

namespace po = boost::program_options;
using namespace std;

const char *input_file_name = "input_data.txt";

struct ParamIsNotSet {
    string param;
    ParamIsNotSet(string param)
        : param(param) {}
};

struct ParamIsNotPositive {
    string param;
    ParamIsNotPositive(string param)
        : param(param) {}
};

void get_double_param(po::variables_map &vm, const char *param, double &val)
{
    if (vm.count(param)) {
        val = vm[param].as<double>();
        if (val < 0)
            throw ParamIsNotPositive(param);
    }
    else
        throw ParamIsNotSet(param);
}

void get_int_param(po::variables_map &vm, const char *param, int &val)
{
    if (vm.count(param)) {
        val = vm[param].as<int>();
        if (val < 0)
            throw ParamIsNotPositive(param);
    }
    else
        throw ParamIsNotSet(param);
}

void get_input_data()
{
    try {

        po::options_description desc("Input data");
        desc.add_options()
            ("L", po::value<double>(), "ll")
            ("T", po::value<double>(), "tt")
            ("a", po::value<double>(), "a")
            ("alpha", po::value<double>(), "alpha")
            ("kappaa", po::value<double>(), "kappaa")
            ("b", po::value<double>(), "b")
            ("beta", po::value<double>(), "beta")
            ("thetab", po::value<double>(), "thetab")
            ("theta0", po::value<double>(), "thetainit")
            ("umin", po::value<double>(), "umin")
            ("umax", po::value<double>(), "umax")
            ("thetad", po::value<double>(), "thetad")
            ("cost_func", po::value<string>(), "cost_func")
            ("init_guess", po::value<string>(), "init_guess")
            ("N", po::value<int>(), "N")
            ("M", po::value<int>(), "tnum")
        ;

        po::variables_map vm;
        ifstream ifs(input_file_name);
        if (!ifs)
        {
            cout << "Can not open input file: " << input_file_name << "\n";
            getch();
            exit(1);
        }
        else
        {
            po::store(parse_config_file(ifs, desc), vm);
            po::notify(vm);
        }

        get_double_param(vm, "L", ll);
        get_double_param(vm, "T", tt);
        get_double_param(vm, "a", a);
        get_double_param(vm, "alpha", alpha);
        get_double_param(vm, "kappaa", kappaa);
        get_double_param(vm, "b", b);
        get_double_param(vm, "beta", beta);
        get_double_param(vm, "thetab", thetab);
        get_double_param(vm, "theta0", thetainit);
        get_double_param(vm, "umin", umin);
        get_double_param(vm, "umax", umax);
        get_double_param(vm, "thetad", thetad);

        get_int_param(vm, "N", N);
        get_int_param(vm, "M", tnum);

        if (vm.count("cost_func")) {
            string s = vm["cost_func"].as<string>();
            if (s == "J1")
                cost_func = COST_FUNC_J1;
            else if (s == "J2")
                cost_func = COST_FUNC_J2;
            else {
                cerr << "cost_func doesn't equal either J1 or J2" << "\n";
                getch();
                exit(1);
            }
        }
        else
            throw ParamIsNotSet("cost_func");

        if (vm.count("init_guess")) {
            string s = vm["init_guess"].as<string>();
            if (s == "umin")
                init_guess_type = INIT_GUESS_UMIN;
            else if (s == "umax")
                init_guess_type = INIT_GUESS_UMAX;
            else {
                cerr << "init_guess doesn't equal either umin or umax" << "\n";
                getch();
                exit(1);
            }
        }
        else
            throw ParamIsNotSet("init_guess");

    }
    catch (ParamIsNotSet e) {
        cerr << "Parameter \"" << e.param << "\" is not specified" << "\n";
        getch();
        exit(1);
    }
    catch (ParamIsNotPositive e) {
        cerr << "Parameter \"" << e.param << "\" is not positive" << "\n";
        getch();
        exit(1);
    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        getch();
        exit(1);
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
        getch();
        exit(1);
    }
}
