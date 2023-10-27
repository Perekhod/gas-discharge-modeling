#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>

// Константы для физических величин
constexpr double PI = 3.141592653589793;
constexpr double K  = 1.38064852e-23;     // постоянная Больцмана
constexpr double E  = 1.60217662e-19;     // элементарный заряд
constexpr double M  = 9.10938356e-31;     // масса электрона

// Класс Particle для представления частицы
class Particle {
private:
    // Поля для хранения координат, скорости, заряда и массы частицы
    double x, y, z;     // координаты в метрах
    double vx, vy, vz;  // скорости в метрах в секунду
    double q;           // заряд в кулонах
    double m;           // масса в килограммах
    bool active;        // поле для отслеживания активности частицы

public:
    // Конструктор класса Particle с параметрами для инициализации полей
    Particle(double x, double y, double z, double vx, double vy, double vz, double q, double m) {
        this->x = x;
        this->y = y;
        this->z = z;
        this->vx = vx;
        this->vy = vy;
        this->vz = vz;
        this->q = q;
        this->m = m;
        this->active = true;
    }

    // Геттеры и сеттеры для доступа к полям класса Particle
    double getX()   { return x; }
    double getY()   { return y; }
    double getZ()   { return z; }
    double getVx()  { return vx;}
    double getVy()  { return vy;}
    double getVz()  { return vz;}
    double getQ()   { return q; }
    double getM()   { return m; }
    bool isActive() { return active;}   // геттер для поля active

    void setX(double x) { this->x = x; }
    void setY(double y) { this->y = y; }
    void setZ(double z) { this->z = z; }
    void setVx(double vx) { this->vx = vx; }
    void setVy(double vy) { this->vy = vy; }
    void setVz(double vz) { this->vz = vz; }
    void setQ(double q) { this->q = q; }
    void setM(double m) { this->m = m; }
    void deactivate()   { active = false;  }  // метод для деактивации частицы
    
    // Метод для обновления координат и скорости частицы в соответствии с электрическим полем и временным шагом
    void update(double Ex, double Ey, double Ez, double dt) {
        // Вычисление ускорения частицы по закону Кулона
        double ax = q * Ex / m;
        double ay = q * Ey / m;
        double az = q * Ez / m;

        // Обновление скорости частицы по закону Ньютона
        vx += ax * dt;
        vy += ay * dt;
        vz += az * dt;

        // Обновление координат частицы по закону движения
        x += vx * dt;
        y += vy * dt;
        z += vz * dt;
    }

    // Метод для определения столкновения частицы с другой частицей или поверхностью
    bool collide(Particle*& other) {
        // Вычисление расстояния между частицами
        double dx = x - other->getX();
        double dy = y - other->getY();
        double dz = z - other->getZ();
        double r = sqrt(dx * dx + dy * dy + dz * dz);

        // Проверка условия столкновения (простая модель сферических частиц)
        return r < 1e-10;
    }

    // Метод для изменения параметров частицы при столкновении с другой частицей или поверхностью
    void bounce(Particle*& other) {
        // Для простоты примера будем считать, что при столкновении частицы обмениваются скоростями и зарядами.

        // Сохранение старых параметров частиц
        double v1x = vx;
        double v1y = vy;
        double v1z = vz;
        double q1 = q;

        // Обмен скоростями и зарядами между частицами
        vx = other->getVx();
        vy = other->getVy();
        vz = other->getVz();
        q = other->getQ();

        // Установка новых скоростей и зарядов для другой частицы
        other->setVx(v1x);
        other->setVy(v1y);
        other->setVz(v1z);
        other->setQ(q1);
    }

    // Метод для определения пересечения траектории частицы с границей моделируемого объема или заданной поверхностью
    bool cross(double Lx, double Ly, double Lz, double Sx, double Sy, double Sz, double Snx, double Sny, double Snz, double dt) {
        // Параметры Lx, Ly, Lz задают размеры моделируемого объема в метрах
        // Параметры Sx, Sy, Sz задают координаты центра заданной поверхности в метрах
        // Параметры Snx, Sny, Snz задают компоненты нормального вектора к заданной поверхности
        // Параметр dt задает временной шаг в секундах

        // Вычисление координат конца траектории частицы после обновления
        double x1 = x + vx * dt;
        double y1 = y + vy * dt;
        double z1 = z + vz * dt;

        // Проверка условия пересечения с границами моделируемого объема
        if (x1 < 0 || x1 > Lx || y1 < 0 || y1 > Ly || z1 < 0 || z1 > Lz) {
            return true;
        }

        // Проверка условия пересечения с заданной поверхностью
        // Для этого используем скалярное произведение векторов
        // Если скалярное произведение вектора нормали к поверхности и вектора траектории частицы отрицательно,
        // то частица пересекает поверхность
        double dx = x - Sx;
        double dy = y - Sy;
        double dz = z - Sz;
        double dx1 = x1 - Sx;
        double dy1 = y1 - Sy;
        double dz1 = z1 - Sz;
        double dot0 = Snx * dx + Sny * dy + Snz * dz;
        double dot1 = Snx * dx1 + Sny * dy1 + Snz * dz1;

        return dot0 * dot1 < 0;
    }
};

// Функция main для моделирования газового разряда
int main() {
    setlocale(LC_ALL, "Ru");

    std::ofstream file("collisions.txt");

    // Задание параметров моделирования
    int N = 1000;           // количество частиц
    double Ex = 1e4;        // напряженность электрического поля по оси x в вольтах на метр
    double Ey = 0;          // напряженность электрического поля по оси y в вольтах на метр
    double Ez = 0;          // напряженность электрического поля по оси z в вольтах на метр
    double dt = 1e-9;       // временной шаг в секундах
    double T = 1e-6;        // длительность моделирования в секундах

    // Задание размеров и формы моделируемого объема
    double Lx = 1e-2;       // длина моделируемого объема по оси x в метрах
    double Ly = 1e-2;       // длина моделируемого объема по оси y в метрах
    double Lz = 1e-2;       // длина моделируемого объема по оси z в метрах

    // Задание характеристик заданной поверхности, с которой создаются новые частицы
    double Sx = 0.5 * Lx;   // координата x центра заданной поверхности в метрах
    double Sy = 0.5 * Ly;   // координата y центра заданной поверхности в метрах
    double Sz = 0;          // координата z центра заданной поверхности в метрах
    double Snx = 0;         // компонента x нормального вектора к заданной поверхности
    double Sny = 0;         // компонента y нормального вектора к заданной поверхности
    double Snz = 1;         // компонента z нормального вектора к заданной поверхности
    double Sd = 1e-3;       // диаметр заданной поверхности в метрах (предполагаем, что она круглая)
    double Sp = 0;          // потенциал заданной поверхности в вольтах

    if (N <= 0) {
        std::cerr << "Ошибка: количество частиц должно быть больше нуля.\n";
        return 1;
    }

    // Создание массива для хранения частиц
    std::vector<Particle*> particles;

    // Создание генератора случайных чисел
    std::random_device rd;
    std::mt19937 gen(rd());

    // Создание равномерного распределения от -1 до 1 для координат и зарядов частиц
    std::uniform_real_distribution<> dis(-1, 1);

    // Создание нормального распределения с нулевым средним и заданной дисперсией для скоростей частиц
    // Дисперсия выбрана так, чтобы соответствовать температуре газа в 300 Кельвинов по формуле (3/2)kT = (1/2)mv^2
    double sigma = sqrt(3 * K * 300 / M);
    std::normal_distribution<> ndis(0, sigma);

    // Цикл для создания и инициализации N частиц на заданной поверхности с случайными координатами и скоростями
    for (int i = 0; i < N; i++) {
        // Генерация случайного угла для определения координат частицы на поверхности
        double theta = 2 * PI * dis(gen);

        // Вычисление координат частицы на поверхности с учетом ее центра и диаметра
        double x = Sx + 0.5 * Sd * cos(theta);
        double y = Sy + 0.5 * Sd * sin(theta);
        double z = Sz;

        // Генерация случайных скоростей частицы в соответствии с распределением Максвелла
        double vx = ndis(gen);
        double vy = ndis(gen);
        double vz = ndis(gen);

        // Вычисление заряда частицы по формуле q = e * (V - Sp) / E, где V - потенциал в точке создания частицы
        // Предполагаем, что потенциал линейно зависит от координаты z: V = Ez * z
        double q = E * (Ez * z - Sp) / E;

        // Создание объекта Particle с заданными параметрами и добавление его в массив
        Particle* p = new Particle(x, y, z, vx, vy, vz, q, M);
        particles.push_back(p);
    }

    int count_part{ 1 }; // счетчик частиц

    // Цикл для моделирования движения частиц в течение времени T с шагом dt
    for (double t = 0; t < T; t += dt) {
        // Цикл для обновления состояния каждой частицы
        for (int i = 0; i < N; i++) {
            // Получение указателя на текущую частицу из массива
            Particle* p = particles[i];

            // Пропуск неактивных частиц
            if (!p->isActive()) continue;  

            // Обновление координат и скорости частицы в соответствии с электрическим полем и временным шагом
            p->update(Ex, Ey, Ez, dt);

            // Проверка условия пересечения траектории частицы с границей моделируемого объема или заданной поверхностью
            if (p->cross(Lx, Ly, Lz, Sx, Sy, Sz, Snx, Sny, Snz, dt)) {
                // Запись координат места пересечения в файл
                file << count_part << " Частица прилипла к поверхности в точке (" << p->getX() << ", " << p->getY() << ", " << p->getZ() << ")\n";
                count_part++;
                // Деактивация частицы при прилипании к поверхности
                p->deactivate();  
            }

            // Цикл для проверки столкновений текущей частицы с другими частицами
            for (int j = i + 1; j < N; j++) {
                // Получение указателя на другую частицу из массива
                Particle* q = particles[j];

                // Проверка условия столкновения между частицами
                if (p->collide(q)) {
                    // Изменение параметров частиц при столкновении
                    p->bounce(q);
                }
            }
        }
    }

    // Объявление переменных для хранения суммарной скорости и энергии частиц
    double vsum{ 0.0 };
    double esum{ 0.0 };

    // Цикл для вычисления суммарной скорости и энергии частиц
    for (int i = 0; i < N; i++) {
        // Получение указателя на текущую частицу из массива
        Particle* p = particles[i];

        // Вычисление модуля скорости частицы
        double v = sqrt(p->getVx() * p->getVx() + p->getVy() * p->getVy() + p->getVz() * p->getVz());

        // Вычисление кинетической энергии частицы
        double e = 0.5 * p->getM() * v * v;

        // Добавление скорости и энергии частицы к сумме
        vsum += v;
        esum += e;
    }

    // Вычисление средней скорости и энергии частиц
    double vavg = vsum / N;
    double eavg = esum / N;

    // Вывод результатов на экран
    std::cout << "Средняя скорость частиц: " << vavg << " м/с\n";
    std::cout << "Средняя энергия частиц: "  << eavg << " Дж\n";

    // Освобождение памяти, занятой частицами
    for (int i = 0; i < N; i++) {
        delete particles[i];
    }

    file.close();
    return 0;
}
