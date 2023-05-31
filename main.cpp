#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <qapplication.h>
#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <QPen>
#include <QColor>

struct Input
{
    int N = 200;
    double T = 0.04;
    double h = T / N;
    double w = 2 * M_PI / T;
    int num = 10;
};

double signal(double x, double T)
{
    double x1, x2, y1, y2;

    if (x >= 0 && x < 0.3 * T)
    {
        x1 = 0; y1 = 5; x2 = 0.3 * T; y2 = -3;
    }
    else if (x >= 0.2 * T && x < 0.5 * T)
    {
        x1 = 0.3 * T; y1 = -3; x2 = 0.5 * T; y2 = -3;
    }
    else if (x >= 0.5 * T && x < 0.8 * T)
    {
        x1 = 0.5 * T; y1 = -3; x2 = 0.8 * T; y2 = 7;
    }
    else if (x >= 0.8 * T && x < T)
    {
        x1 = 0.8 * T; y1 = 7; x2 = T; y2 = 0;
    }
    else
    {
        return 0;
    }
    return (x - x1) * (y2 - y1) / (x2 - x1) + y1;
}

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    QwtPlot plot;
    plot.setTitle("Signal");
    plot.setCanvasBackground(Qt::white);

    QwtPlotCurve signalCurve("Signal");
    signalCurve.setRenderHint(QwtPlotItem::RenderAntialiased);
    QColor waveColor(0, 0, 255); // Синій колір
    QPen pen(waveColor);
    pen.setWidth(2);
    pen.setStyle(Qt::SolidLine);
    pen.setCosmetic(true);
    pen.setCapStyle(Qt::FlatCap);
    pen.setJoinStyle(Qt::RoundJoin);
    signalCurve.setPen(pen);

    Input input;

    std::vector<double> t(input.N), s(input.N);
    for (int i = 0; i < input.N; ++i)
    {
        t[i] = i * input.h;
        s[i] = signal(t[i], input.T) + 0.2 * sin(20 * t[i]);
    }
    signalCurve.setSamples(t.data(), s.data(), input.N);
    signalCurve.attach(&plot);

    plot.replot();
    plot.show();

    // Part III
    std::ofstream coefficientsFile("coefficients.txt");
    if (coefficientsFile.is_open()) {
        for (int n = 0; n <= input.num; ++n) {
            double realSum = 0.0, imagSum = 0.0;
            for (int i = 0; i < input.N; ++i) {
                realSum += s[i] * cos(n * input.w * t[i]);
                imagSum += s[i] * sin(n * input.w * t[i]);
            }
            double coefficient = (1.0 / input.N) * (realSum * cos(n * input.w * t[0]) + imagSum * sin(n * input.w * t[0]));
            coefficientsFile << coefficient << "\n";
        }
        coefficientsFile.close();
    }
    else {
        std::cout << "Error opening coefficients file." << std::endl;
        return 1;
    }

    // Part IV
    std::ofstream sumFile("sum.txt");
    if (sumFile.is_open()) {
        for (int i = 0; i < input.N; ++i) {
            double sum = 0.0;
            for (int n = 0; n <= input.num; ++n) {
                double coefficient;
                std::ifstream coefficientsFile("coefficients.txt");
                if (coefficientsFile.is_open()) {
                    for (int j = 0; j <= n; ++j) {
                        coefficientsFile >> coefficient;
                    }
                    coefficientsFile.close();
                    sum += coefficient * cos(n * input.w * t[i]);
                }
                else {
                    std::cout << "Error opening coefficients file." << std::endl;
                    return 1;
                }
            }
            sumFile << sum << "\n";
        }
        sumFile.close();
    }
    else {
        std::cout << "Error opening sum file." << std::endl;
        return 1;
    }

    // Part V:
    QwtPlot sumPlot;
    sumPlot.setTitle("Fourier Series Sum");
    sumPlot.setCanvasBackground(Qt::white);

    QwtPlotCurve sumCurve("Sum");
    sumCurve.setRenderHint(QwtPlotItem::RenderAntialiased);
    sumCurve.setPen(QPen(Qt::red));
    std::vector<double> sumValues(input.N);
    std::ifstream sumFileRead("sum.txt");
    if (sumFileRead.is_open()) {
        for (int i = 0; i < input.N; ++i) {
            sumFileRead >> sumValues[i];
        }
        sumFileRead.close();
    }
    else {
        std::cout << "Error opening sum file for reading." << std::endl;
        return 1;
    }

    sumCurve.setSamples(t.data(), sumValues.data(), input.N);
    sumCurve.attach(&sumPlot);

    sumPlot.replot();
    sumPlot.show();

    return app.exec();
}
