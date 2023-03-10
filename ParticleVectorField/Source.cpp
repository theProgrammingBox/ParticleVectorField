#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::nanoseconds;
using std::chrono::microseconds;

/*
IMPORTANT LESSONS
1. With euler method, a substep of 10 is needed to be stable
2. With RK4 method, a substep of 920 is needed to be stable
3. For about 4x the computation, RK4 allows you to have about 92x the stability
4. However, you need more storage for RK4 than euler method
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
	static constexpr int numParticles = 100000;
	double* x = new double[numParticles];
	double* y = new double[numParticles];
	double* z = new double[numParticles];
	double orginX, orginY;
	
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

		orginX = ScreenWidth() * 0.5;
		orginY = ScreenHeight() * 0.5;
		
		return true;
	}

	int particlesInScreen()
	{
		int count = 0;
		for (int i = numParticles; i--;)
			count += x[i] * 2 + orginX > 0 && x[i] * 2 + orginX < ScreenWidth() && y[i] * 2 + orginY > 0 && y[i] * 2 + orginY < ScreenHeight();
		return count;
	}

	void Euler(double dt)
	{
		double dx, dy, dz;

		for (int i = numParticles; i--;)
		{
			halvorsenAttractor(x[i], y[i], z[i], dx, dy, dz);

			x[i] += dx * dt;
			y[i] += dy * dt;
			z[i] += dz * dt;
		}
	}

	void RungeKutta(double dt)
	{
		double tempX, tempY, tempZ, dx, dy, dz, savedX, savedY, savedZ;
		double applied[4] = { 0 * dt, 0.5 * dt, 0.5 * dt, 1 * dt };
		double summed[4] = { 0.166666666667 * dt, 0.333333333333 * dt, 0.333333333333 * dt, 0.166666666667 * dt };

		for (uint32_t i = numParticles; i--;)
		{
			savedX = x[i];
			savedY = y[i];
			savedZ = z[i];
			dx = dy = dz = 0;
			
			for (uint32_t j = 0; j < 4; j++)
			{
				tempX = savedX + applied[j] * dx;
				tempY = savedY + applied[j] * dy;
				tempZ = savedZ + applied[j] * dz;
				halvorsenAttractor(tempX, tempY, tempZ, dx, dy, dz);
				x[i] += summed[j] * dx;
				y[i] += summed[j] * dy;
				z[i] += summed[j] * dz;
			}
		}
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		double dt = 0.1;
		int subSteps;
		bool rungeKutta = true;
		
		if (rungeKutta)
		{
			subSteps = 1;
			dt /= subSteps;
			for (int i = subSteps; i--;)
				RungeKutta(dt);
		}
		else
		{
			subSteps = 90;
			dt /= subSteps;
			for (int i = subSteps; i--;)
				Euler(dt);
		}
		
		Clear(olc::BLACK);

		for (int i = numParticles; i--;)
			Draw(x[i] * 20 + orginX, y[i] * 20 + orginY);
		
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