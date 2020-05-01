/**
 *
 * "A Collection of Useful C++ Classes for Digital Signal Processing"
 * By Vinnie Falco and Bernd Porr
 *
 * Official project location:
 * https://github.com/berndporr/iir1
 *
 * See Documentation.cpp for contact information, notes, and bibliography.
 * 
 * -----------------------------------------------------------------
 *
 * License: MIT License (http://www.opensource.org/licenses/mit-license.php)
 * Copyright (c) 2009 by Vinnie Falco
 * Copyright (c) 2011 by Bernd Porr
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 **/

#ifndef IIR1_RBJ_H
#define IIR1_RBJ_H

#include "Common.h"
#include "Biquad.h"
#include "State.h"

namespace Iir
{

/**
 * Filter realizations based on Robert Bristol-Johnson formulae:
 *
 * http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
 *
 * These are all 2nd order filters which are tuned with the Q (or Quality factor).
 * The Q factor causes a resonance at the cutoff frequency. The higher the Q
 * factor the higher the responance. If 0.5 < Q < 1/sqrt(2) then there is no resonance peak.
 * Above 1/sqrt(2) the peak becomes more and more pronounced. For bandpass and stopband
 * the Q factor is replaced by the width of the filter. The higher Q the more narrow
 * the bandwidth of the notch or bandpass.
 *
 **/

#define ONESQRT2 (1/sqrt(2))

	namespace RBJ
	{

		/** 
			 * The base class of all RBJ filters
			 **/
		template<typename Sample=double>
		struct DllExport RBJbase : Biquad
		{
		public:
			/// filter operation
			inline Sample filter(const Sample& s)
			{
				return (state.filter(s, *this));
			}
			/// resets the delay lines to zero
			void reset()
			{
				state.reset();
			}
			/// gets the delay lines (=state) of the filter
			const DirectFormI <Sample>& getState()
			{
				return state;
			}
		private:
			DirectFormI <Sample> state;
		};

		/**
			 * Lowpass.
			 **/
		template<typename Sample=double>
		struct DllExport LowPass : RBJbase<Sample>
		{
			/**
					 * Calculates the coefficients
					 * \param sampleRate Sampling rate
					 * \param cutoffFrequency Cutoff frequency
					 * \param q Q factor determines the resonance peak at the cutoff.
					 **/
			void setup(double sampleRate,
				double cutoffFrequency,
				double q = ONESQRT2)
			{
				double w0 = 2 * doublePi * cutoffFrequency / sampleRate;
				double cs = cos(w0);
				double sn = sin(w0);
				double AL = sn / (2 * q);
				double b0 = (1 - cs) / 2;
				double b1 = 1 - cs;
				double b2 = (1 - cs) / 2;
				double a0 = 1 + AL;
				double a1 = -2 * cs;
				double a2 = 1 - AL;
				Biquad::setCoefficients(a0, a1, a2, b0, b1, b2);
			}
		};

		/**
			 * Highpass.
			 **/
		template<typename Sample=double>
		struct DllExport HighPass : RBJbase<Sample>
		{
			/**
					 * Calculates the coefficients
					 * \param sampleRate Sampling rate
					 * \param cutoffFrequency Cutoff frequency
					 * \param q Q factor determines the resonance peak at the cutoff.
					 **/
			void setup(double sampleRate,
				double cutoffFrequency,
				double q = ONESQRT2)
			{
				double w0 = 2 * doublePi * cutoffFrequency / sampleRate;
				double cs = cos(w0);
				double sn = sin(w0);
				double AL = sn / (2 * q);
				double b0 = (1 + cs) / 2;
				double b1 = -(1 + cs);
				double b2 = (1 + cs) / 2;
				double a0 = 1 + AL;
				double a1 = -2 * cs;
				double a2 = 1 - AL;
				Biquad::setCoefficients(a0, a1, a2, b0, b1, b2);
			}
		};

		/**
			 * Bandpass with constant skirt gain
			 **/
		template<typename Sample=double>
		struct DllExport BandPass1 : RBJbase<Sample>
		{
			/**
					 * Calculates the coefficients
					 * \param sampleRate Sampling rate
					 * \param centerFrequency Center frequency of the bandpass
					 * \param bandWidth Bandwidth in octaves
					 **/
			void setup(double sampleRate,
				double centerFrequency,
				double bandWidth)
			{
				double w0 = 2 * doublePi * centerFrequency / sampleRate;
				double cs = cos(w0);
				double sn = sin(w0);
				double AL = sn / (2 * bandWidth);
				double b0 = bandWidth * AL;// sn / 2;
				double b1 = 0;
				double b2 = -bandWidth * AL;//-sn / 2;
				double a0 = 1 + AL;
				double a1 = -2 * cs;
				double a2 = 1 - AL;
				Biquad::setCoefficients(a0, a1, a2, b0, b1, b2);
			}
		};

		/**
			 * Bandpass with constant 0 dB peak gain
			 **/
		template<typename Sample=double>
		struct DllExport BandPass2 : RBJbase<Sample>
		{
			/**
					 * Calculates the coefficients
					 * \param sampleRate Sampling rate
					 * \param centerFrequency Center frequency of the bandpass
					 * \param bandWidth Bandwidth in octaves
					 **/
			void setup(double sampleRate,
				double centerFrequency,
				double bandWidth)
			{
				double w0 = 2 * doublePi * centerFrequency / sampleRate;
				double cs = cos(w0);
				double sn = sin(w0);
				double AL = sn / (2 * bandWidth);
				double b0 = AL;
				double b1 = 0;
				double b2 = -AL;
				double a0 = 1 + AL;
				double a1 = -2 * cs;
				double a2 = 1 - AL;
				Biquad::setCoefficients(a0, a1, a2, b0, b1, b2);
			}
		};

		/**
			 * Bandstop filter. Warning: the bandwidth might not be accurate
			 * for narrow notches.
			 **/
		template<typename Sample=double>
		struct DllExport BandStop : RBJbase<Sample>
		{
			/**
					 * Calculates the coefficients
					 * \param sampleRate Sampling rate
					 * \param centerFrequency Center frequency of the bandstop
					 * \param bandWidth Bandwidth in octaves
					 **/
			void setup(double sampleRate,
				double centerFrequency,
				double bandWidth)
			{
				double w0 = 2 * doublePi * centerFrequency / sampleRate;
				double cs = cos(w0);
				double sn = sin(w0);
				double AL = sn / (2 * bandWidth);
				double b0 = 1;
				double b1 = -2 * cs;
				double b2 = 1;
				double a0 = 1 + AL;
				double a1 = -2 * cs;
				double a2 = 1 - AL;
				Biquad::setCoefficients(a0, a1, a2, b0, b1, b2);
			}
		};

		/**
			 * Bandstop with Q factor: the higher the Q factor the more narrow is
			 * the notch.
			 * However, a narrow notch has a long impulse response ( = ringing)
			 * and numerical problems might prevent perfect damping. Practical values
			 * of the Q factor are about Q = 10 to 20. In terms of the design
			 * the Q factor defines the radius of the
			 * poles as r = exp(- pi*(centerFrequency/sampleRate)/q_factor) whereas
			 * the angles of the poles/zeros define the bandstop frequency. The higher
			 * Q the closer r moves towards the unit circle.
			 **/
		template<typename Sample=double>
		struct DllExport IIRNotch : RBJbase<Sample>
		{
			/**
					 * Calculates the coefficients
					 * \param sampleRate Sampling rate
					 * \param centerFrequency Center frequency of the notch
					 * \param q_factor Q factor of the notch (1 to ~20)
					 **/
			void setup(double sampleRate,
				double centerFrequency,
				double q_factor = 10)
			{

				double w0 = 2 * doublePi * centerFrequency / sampleRate;
				double cs = cos(w0);
				double r = exp(-(w0 / 2) / q_factor);
				double b0 = 1;
				double b1 = -2 * cs;
				double b2 = 1;
				double a0 = 1;
				double a1 = -2 * r * cs;
				double a2 = r * r;
				Biquad::setCoefficients(a0, a1, a2, b0, b1, b2);
			}
		};

		/**
			 * Low shelf: 0db in the stopband and gainDb in the passband.
			 **/
		template<typename Sample=double>
		struct DllExport LowShelf : RBJbase<Sample>
		{
			/**
			 * Calculates the coefficients
			 * \param sampleRate Sampling rate
			 * \param cutoffFrequency Cutoff frequency
			 * \param gainDb Gain in the passband
					 * \param shelfSlope Slope between stop/passband. 1 = as steep as it can.
			 **/
			void setup(double sampleRate,
				double cutoffFrequency,
				double gainDb,
				double shelfSlope = 1)
			{
				double A = pow(10, gainDb / 40);
				double w0 = 2 * doublePi * cutoffFrequency / sampleRate;
				double cs = cos(w0);
				double sn = sin(w0);
				double AL = sn / 2 * ::std::sqrt((A + 1 / A) * (1 / shelfSlope - 1) + 2);
				double sq = 2 * sqrt(A) * AL;
				double b0 = A * ((A + 1) - (A - 1) * cs + sq);
				double b1 = 2 * A * ((A - 1) - (A + 1) * cs);
				double b2 = A * ((A + 1) - (A - 1) * cs - sq);
				double a0 = (A + 1) + (A - 1) * cs + sq;
				double a1 = -2 * ((A - 1) + (A + 1) * cs);
				double a2 = (A + 1) + (A - 1) * cs - sq;
				Biquad::setCoefficients(a0, a1, a2, b0, b1, b2);
			}
		};

		/**
			 * High shelf: 0db in the stopband and gainDb in the passband.
			 **/
		template<typename Sample=double>
		struct DllExport HighShelf : RBJbase<Sample>
		{
			/**
			 * Calculates the coefficients
			 * \param sampleRate Sampling rate
			 * \param cutoffFrequency Cutoff frequency
			 * \param gainDb Gain in the passband
					 * \param shelfSlope Slope between stop/passband. 1 = as steep as it can.
			 **/
			void setup(double sampleRate,
				double cutoffFrequency,
				double gainDb,
				double shelfSlope = 1)
			{
				double A = pow(10, gainDb / 40);
				double w0 = 2 * doublePi * cutoffFrequency / sampleRate;
				double cs = cos(w0);
				double sn = sin(w0);
				double AL = sn / 2 * ::std::sqrt((A + 1 / A) * (1 / shelfSlope - 1) + 2);
				double sq = 2 * sqrt(A) * AL;
				double b0 = A * ((A + 1) + (A - 1) * cs + sq);
				double b1 = -2 * A * ((A - 1) + (A + 1) * cs);
				double b2 = A * ((A + 1) + (A - 1) * cs - sq);
				double a0 = (A + 1) - (A - 1) * cs + sq;
				double a1 = 2 * ((A - 1) - (A + 1) * cs);
				double a2 = (A + 1) - (A - 1) * cs - sq;
				Biquad::setCoefficients(a0, a1, a2, b0, b1, b2);
			}
		};

		/**
			 * Band shelf: 0db in the stopband and gainDb in the passband.
			 **/
		template<typename Sample=double>
		struct DllExport BandShelf : RBJbase<Sample>
		{
			/**
			 * Calculates the coefficients
			 * \param sampleRate Sampling rate
			 * \param centerFrequency frequency
			 * \param gainDb Gain in the passband
					 * \param bandWidth Bandwidth in octaves
			 **/
			void setup(double sampleRate,
				double centerFrequency,
				double gainDb,
				double bandWidth)
			{
				double A = pow(10, gainDb / 40);
				double w0 = 2 * doublePi * centerFrequency / sampleRate;
				double cs = cos(w0);
				double sn = sin(w0);
				double AL = sn * sinh(doubleLn2 / 2 * bandWidth * w0 / sn);
				if (Iir::is_nan(AL))
					throw std::invalid_argument("No solution available for these parameters.\n");
				double b0 = 1 + AL * A;
				double b1 = -2 * cs;
				double b2 = 1 - AL * A;
				double a0 = 1 + AL / A;
				double a1 = -2 * cs;
				double a2 = 1 - AL / A;
				Biquad::setCoefficients(a0, a1, a2, b0, b1, b2);
			}
		};

		template<typename Sample=double>
		struct DllExport AllPass : RBJbase<Sample>
		{
			void setup(double sampleRate,
				double phaseFrequency,
				double q)
			{
				double w0 = 2 * doublePi * phaseFrequency / sampleRate;
				double cs = cos(w0);
				double sn = sin(w0);
				double AL = sn / (2 * q);
				double b0 = 1 - AL;
				double b1 = -2 * cs;
				double b2 = 1 + AL;
				double a0 = 1 + AL;
				double a1 = -2 * cs;
				double a2 = 1 - AL;
				Biquad::setCoefficients(a0, a1, a2, b0, b1, b2);
			}

		};

	}

}

#endif
