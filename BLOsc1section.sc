BLOsc1section : PureUGen  {
	*ar {
		arg
		freq = 440.0,

		loHarmonics1 = 1.0,
		hiHarmonics1 = 5.0,
		slope1 = 1.0,
		evenOddRatio1 = 1.0,
		spread1 = 1,
		spreadCompensation = 1,

		mul = 1.0,
		add = 0.0;

		^this.multiNew('audio',
		freq,

		loHarmonics1,
		hiHarmonics1,
		slope1,
		evenOddRatio1,
		spread1,
		spreadCompensation,
		).madd(mul, add)
	}

	*kr { arg
		freq = 1.0,

		loHarmonics1 = 1.0,
		hiHarmonics1 = 5.0,
		slope1 = 1.0,
		evenOddRatio1 = 1.0,
		spread1 = 1,
		spreadCompensation = 1,

		mul = 1.0,
		add = 0.0;

		^this.multiNew('control',
		freq,

		loHarmonics1,
		hiHarmonics1,
		slope1,
		evenOddRatio1,
		spread1,
		spreadCompensation,
		).madd(mul, add)
	}
}
