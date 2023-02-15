#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::nanoseconds;
using std::chrono::microseconds;

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

// Override base class with your custom functionality
class Example : public olc::PixelGameEngine
{
public:
	static constexpr int numParticles = 10000;
	double x[numParticles];
	double y[numParticles];
	double z[numParticles];

	double orginX, orginY;
	
	Example()
	{
		sAppName = "Example";
	}
	
	void halvorsenAttractor(double x, double y, double z, double& dx, double& dy, double& dz)
	{
		const double a = 1.4;

		x = x * 0.1;
		y = y * 0.1;
		z = z * 0.1;

		dx = -a * x - 4 * y - 4 * z - y * y;
		dy = -a * y - 4 * z - 4 * x - z * z;
		dz = -a * z - 4 * x - 4 * y - x * x;
	}

	
	bool OnUserCreate() override
	{
		for (int i = 0; i < numParticles; i++) {
			x[i] = GLOBAL::random.Rfloat(-80, 20);
			y[i] = GLOBAL::random.Rfloat(-80, 20);
			z[i] = GLOBAL::random.Rfloat(-80, 20);
		}

		orginX = ScreenWidth() / 2;
		orginY = ScreenHeight() / 2;
		
		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override
	{
		double dt = 0.1f;
		double dx1[numParticles];
		double dy1[numParticles];
		double dz1[numParticles];
		double dx2[numParticles];
		double dy2[numParticles];
		double dz2[numParticles];
		double dx3[numParticles];
		double dy3[numParticles];
		double dz3[numParticles];
		double dx4[numParticles];
		double dy4[numParticles];
		double dz4[numParticles];
		
		// Runge-Kutta (RK4)
		for (int i = numParticles; i--;)
		{
			halvorsenAttractor(x[i], y[i], z[i], dx1[i], dy1[i], dz1[i]);
			halvorsenAttractor(x[i] + dx1[i] * dt / 2, y[i] + dy1[i] * dt / 2, z[i] + dz1[i] * dt / 2, dx2[i], dy2[i], dz2[i]);
			halvorsenAttractor(x[i] + dx2[i] * dt / 2, y[i] + dy2[i] * dt / 2, z[i] + dz2[i] * dt / 2, dx3[i], dy3[i], dz3[i]);
			halvorsenAttractor(x[i] + dx3[i] * dt, y[i] + dy3[i] * dt, z[i] + dz3[i] * dt, dx4[i], dy4[i], dz4[i]);
			
			x[i] += (dx1[i] + 2 * dx2[i] + 2 * dx3[i] + dx4[i]) * dt / 6;
			y[i] += (dy1[i] + 2 * dy2[i] + 2 * dy3[i] + dy4[i]) * dt / 6;
			z[i] += (dz1[i] + 2 * dz2[i] + 2 * dz3[i] + dz4[i]) * dt / 6;
		}
		
		Clear(olc::BLACK);
		for (int i = numParticles; i--;)
			Draw(x[i] * 2 + orginX, y[i] * 2 + orginY);
		
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