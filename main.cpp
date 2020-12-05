#include "tiny_pga.h"
#include <iostream>

using namespace tiny_pga;

int main()
{
  tiny_pga::Plane plane{1.f, 2.f, 3.f, 4.f};

  tiny_pga::Point point{{}, {}, {}, 1.f, 2.f, 3.f, 0.F};

  auto a = plane * point;

  return 0;
}