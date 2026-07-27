#include <kids/core/kidsdata.h>
#include <kids/toltec/timestream.h>

int main()
{
    auto *reader = &kids::toltec::read_raw_timestream_slice;
    return reader == nullptr;
}
