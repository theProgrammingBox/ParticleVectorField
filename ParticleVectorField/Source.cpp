#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::nanoseconds;
using std::chrono::microseconds;

/*
IMPORTANT LESSONS
1. With euler method, a time step less than 0.0016 is needed to be stable
2. With RK4 method, a time step less than 0.1 is needed to be stable
3. For about 4x the computation, RK4 allows you to have about 60x the time step
4. You need more storage for RK4 than euler method
*/

class Random
{
public:
	Random(uint32_t seed = 0)	// seed the random number generator
	{
		state[0] = Hash((uint8_t*)&seed, 4, seed);
		state[1] = Hash((uint8_t*)&seed, 4, state[0]);
	}

	static uint32_t MakeSeed(uint32_t seed = 0)	// make seed from time and seed
	{
		uint32_t result = seed;
		result = Hash((uint8_t*)&result, 4, nanosecond());
		result = Hash((uint8_t*)&result, 4, microsecond());
		return result;
	}

	void Seed(uint32_t seed = 0)	// seed the random number generator
	{
		state[0] = Hash((uint8_t*)&seed, 4, seed);
		state[1] = Hash((uint8_t*)&seed, 4, state[0]);
	}

	uint32_t Ruint32()	// XORSHIFT128+
	{
		uint64_t a = state[0];
		uint64_t b = state[1];
		state[0] = b;
		a ^= a << 23;
		state[1] = a ^ b ^ (a >> 18) ^ (b >> 5);
		return uint32_t((state[1] + b) >> 16);
	}

	float Rfloat(float min = 0, float max = 1) { return min + (max - min) * Ruint32() * 2.3283064371e-10; }

	static uint32_t Hash(const uint8_t* key, size_t len, uint32_t seed = 0)	// MurmurHash3
	{
		uint32_t h = seed;
		uint32_t k;
		for (size_t i = len >> 2; i; i--) {
			memcpy(&k, key, 4);
			key += 4;
			h ^= murmur_32_scramble(k);
			h = (h << 13) | (h >> 19);
			h = h * 5 + 0xe6546b64;
		}
		k = 0;
		for (size_t i = len & 3; i; i--) {
			k <<= 8;
			k |= key[i - 1];
		}
		h ^= murmur_32_scramble(k);
		h ^= len;
		h ^= h >> 16;
		h *= 0x85ebca6b;
		h ^= h >> 13;
		h *= 0xc2b2ae35;
		h ^= h >> 16;
		return h;
	}

private:
	uint64_t state[2];

	static uint32_t murmur_32_scramble(uint32_t k) {
		k *= 0xcc9e2d51;
		k = (k << 15) | (k >> 17);
		k *= 0x1b873593;
		return k;
	}

	static uint32_t nanosecond() { return duration_cast<nanoseconds>(high_resolution_clock::now().time_since_epoch()).count(); }
	static uint32_t microsecond() { return duration_cast<microseconds>(high_resolution_clock::now().time_since_epoch()).count(); }
};

namespace GLOBAL
{
	Random random(Random::MakeSeed(0));
}

class Example : public olc::PixelGameEngine
{
public:
	static constexpr int numParticles = 10000;
	double x[numParticles];
	double y[numParticles];
	double z[numParticles];
	double orginX, orginY;
	double minSpeed, maxSpeed;
	
	Example()
	{
		sAppName = "Example";
	}
	
	void halvorsenAttractor(double x, double y, double z, double& dx, double& dy, double& dz)
	{
		dx = -1.4 * x - 4 * y - 4 * z - y * y;
		dy = -1.4 * y - 4 * z - 4 * x - z * z;
		dz = -1.4 * z - 4 * x - 4 * y - x * x;
	}

	
	bool OnUserCreate() override
	{
		for (int i = 0; i < numParticles; i++)
		{
			x[i] = GLOBAL::random.Rfloat(-8, 2);
			y[i] = GLOBAL::random.Rfloat(-8, 2);
			z[i] = GLOBAL::random.Rfloat(-8, 2);
		}

		orginX = ScreenWidth() / 2;
		orginY = ScreenHeight() / 2;
		
		return true;
	}

	int particlesInScreen()
	{
		int count = 0;
		for (int i = numParticles; i--;)
			count += x[i] * 2 + orginX > 0 && x[i] * 2 + orginX < ScreenWidth() && y[i] * 2 + orginY > 0 && y[i] * 2 + orginY < ScreenHeight();
		return count;
	}

	olc::Pixel mapToRainbow(double d)	// 0 - 1
	{
		d *= 6;
		double r = (d > 4) ? std::max(0.0, std::min(1.0, 6 - d)) : std::max(0.0, std::min(1.0, d - 2));
		double g = (d > 2) ? std::max(0.0, std::min(1.0, 4 - d)) : std::max(0.0, std::min(1.0, d));
		double b = (d > 3) ? std::max(0.0, std::min(1.0, d - 4)) : std::max(0.0, std::min(1.0, 2 - d));

		return olc::Pixel(r * 0xff, g * 0xff, b * 0xff);
	}

	olc::Pixel mapToBAndW(double d) { return olc::Pixel(d * 0xff, d * 0xff, d * 0xff); }

	void RungeKutta(double dt)
	{
		Clear(olc::BLACK);
		
		double dx1, dy1, dz1;
		double dx2, dy2, dz2;
		double dx3, dy3, dz3;
		double dx4, dy4, dz4;
		double range = maxSpeed - minSpeed;
		double speed;
		double tempMinSpeed = DBL_MAX;
		double tempMaxSpeed = DBL_MIN;
		
		for (int i = numParticles; i--;)
		{
			halvorsenAttractor(x[i], y[i], z[i], dx1, dy1, dz1);
			halvorsenAttractor(x[i] + dx1 * dt * 0.5, y[i] + dy1 * dt * 0.5, z[i] + dz1 * dt * 0.5, dx2, dy2, dz2);
			halvorsenAttractor(x[i] + dx2 * dt * 0.5, y[i] + dy2 * dt * 0.5, z[i] + dz2 * dt * 0.5, dx3, dy3, dz3);
			halvorsenAttractor(x[i] + dx3 * dt, y[i] + dy3 * dt, z[i] + dz3 * dt, dx4, dy4, dz4);

			dx1 = (dx1 + 2 * dx2 + 2 * dx3 + dx4) * dt * 0.166666666667;
			dy1 = (dy1 + 2 * dy2 + 2 * dy3 + dy4) * dt * 0.166666666667;
			dz1 = (dz1 + 2 * dz2 + 2 * dz3 + dz4) * dt * 0.166666666667;
			
			x[i] += dx1;
			y[i] += dy1;
			z[i] += dz1;

			speed = std::sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
			tempMinSpeed = std::min(tempMinSpeed, speed);
			tempMaxSpeed = std::max(tempMaxSpeed, speed);
			
			//Draw(x[i] * 20 + orginX, y[i] * 20 + orginY, mapToRainbow((speed - minSpeed) / range));
			Draw(x[i] * 20 + orginX, y[i] * 20 + orginY, mapToBAndW((speed - minSpeed) / range));
		}
		
		minSpeed = tempMinSpeed;
		maxSpeed = tempMaxSpeed;
	}

	void Euler(double dt)
	{
		Clear(olc::BLACK);
		
		double dx, dy, dz;
		double range = maxSpeed - minSpeed;
		double speed;
		double tempMinSpeed = DBL_MAX;
		double tempMaxSpeed = DBL_MIN;

		for (int i = numParticles; i--;)
		{
			halvorsenAttractor(x[i], y[i], z[i], dx, dy, dz);

			x[i] += dx * dt;
			y[i] += dy * dt;
			z[i] += dz * dt;

			speed = std::sqrt(dx * dx + dy * dy + dz * dz);
			tempMinSpeed = std::min(tempMinSpeed, speed);
			tempMaxSpeed = std::max(tempMaxSpeed, speed);
			
			//Draw(x[i] * 20 + orginX, y[i] * 20 + orginY, mapToRainbow((speed - minSpeed) / range));
			Draw(x[i] * 20 + orginX, y[i] * 20 + orginY, mapToBAndW((speed - minSpeed) / range));
		}
		
		minSpeed = tempMinSpeed;
		maxSpeed = tempMaxSpeed;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		double dt = 0.1;
		bool rungeKutta = true;
		
		if (rungeKutta)
		{
			RungeKutta(dt);
		}
		else
		{
			Euler(dt);
		}
		
		DrawString(10, 10, "Particles in screen: " + std::to_string(particlesInScreen()) + " / " + std::to_string(numParticles), olc::WHITE);
		
		return true;
	}
};

int main()
{
	Example demo;
	if (demo.Construct(1280, 720, 1, 1))
		demo.Start();
	return 0;
}