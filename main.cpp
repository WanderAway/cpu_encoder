#include "encoder.h"
#include <chrono>

int main(int argc, char* argv[])
{
  if (argc < 3)
  {
    cout << "Usage: " << argv[0] << " <parallization mode> <qp>\n";
    cout << "       mode 0 - no parallel processing\n";
    cout << "       mode 1 - block level parallization\n";
    cout << "       mode 2 - frame level parallization\n";
    return 0;
  }
  u8 pm = atoi(argv[1]);
  u8 qp = atoi(argv[2]);
  encoder e(1, vector2d(352,288),16,4,pm,qp);

  ifstream infile("foreman.y");
  ofstream outfile("generated.y");

  auto start = chrono::system_clock::now();

  e.encode(infile, outfile);

  auto end = chrono::system_clock::now();
  auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);

  cout <<"Runtime: " << elapsed.count() << "ms" << endl;
  
  infile.close();
  outfile.close();

  return 0;
}
