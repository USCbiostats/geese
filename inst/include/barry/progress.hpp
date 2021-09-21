#ifndef BARRY_PROGRESS_HPP
#define BARRY_PROGRESS_HPP

#ifndef BARRY_PROGRESS_BAR_WIDTH
#define BARRY_PROGRESS_BAR_WIDTH 80
#endif

class Progress {
private:
    int width;
    int n;
    int step;
    int n_steps;
    int step_size;
    int n_bars;
public:

    Progress(int n_, int width_);
    ~Progress() {};

    void next();
    void end();

};

inline Progress::Progress(int n_, int width_) {


    width   = std::max(1, width_ - 7);
    n       = n_;
    step    = 0;
    n_steps = (n > width) ?
        static_cast<int>(floor(static_cast<double>(n) / width)) : n;

    step_size  = static_cast<int>(n / n_steps);

    n_bars = static_cast<int>(std::max(1.0, floor(width / static_cast<double>(n_steps))));

}

inline void Progress::next() {

    if (!(step++ % step_size))
    {

        if (width <= 1)
            return;

        for (int j = 0; j < n_bars; ++j)
            printf_barry("|");
    }

}

inline void Progress::end() {

    int reminder = static_cast<int>(width) - n_bars * n_steps;

    if ((width > 1) & (reminder > 0))
    {

        for (int j = 0; j < reminder; ++j)
            printf_barry("|");

    }
    
    printf_barry(" done.\n");

}

#endif