
#include <savvy/reader.hpp>

#include <iostream>
#include <cstdlib>

int main(int argc, char** argv)
{

  if (argc == 2)
  {
    savvy::reader rdr(argv[1], savvy::fmt::gt);

    if (rdr)
    {
      for (auto it = rdr.samples().begin(); it != rdr.samples().end(); ++it)
        std::cout << (*it) << std::endl;
      return EXIT_SUCCESS;
    }
  }

  return EXIT_FAILURE;
}