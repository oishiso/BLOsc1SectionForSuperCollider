// SynthDef
(
SynthDef(\BLOsc1sectionExample, {
	arg
	freq = 220.0,

	loHarmonics1 = 1.0,
	hiHarmonics1 = 6.0,
	slope1 = 1.0,
	evenOddRatio1 = 1.0,
	spread1 = 1,
	spreadCompensation = 1,

	amp = 1.0
	;

	Out.ar(0,
		BLOsc1section.ar(
			freq: freq,

			loHarmonics1: loHarmonics1,
			hiHarmonics1: hiHarmonics1,
			slope1: slope1,
			evenOddRatio1: evenOddRatio1,
			spread1: spread1,

			spreadCompensation: spreadCompensation,

			mul: amp))
}).add;
)

// Spectrum Analyzer
(
~openAnalyzer = {
	var window, func;
	window = Window("Analyzer", Rect(100, 100, 511, 300)); // width should be 511
	func = FreqScopeView(window, window.view.bounds);
	func.active_(true); // turn it on the first time;
	func.freqMode_(1); //log scale
	window.onClose_({ func.kill }); // you must have this
	window.front;
}.value
)

FreqScope.new

(
{
	BLOsc1section.ar(
		freq: SinOsc.ar(0.01).range(100,200),

		loHarmonics1: 1.0,
		hiHarmonics1: LFNoise2.ar(0.4).range(3.0, 5.0),
		slope1: LFNoise2.ar(0.3).range(0.1, 1.9),
		evenOddRatio1: LFNoise2.ar(0.2).range(0.0, 1),
		spread1: 2,

		spreadCompensation: 1,

		mul: 1.0,
		add: 0.0
	)
}.play
)


// Demonstration with Gui
(
var window, sound,

knob_freq,

knob_loHarmonics1,
knob_hiHarmonics1,
knob_slope1,
knob_evenOddRatio1,
knob_spread1,


knob_spreadCompensation,
knob_amp;


window = Window("Band Limited Oscillator", Rect(left: 100, top: 100, width: 650, height: 420));
window.view.decorator = FlowLayout(window.view.bounds);
StaticText(window,500 @ 40).string_("Band Limited Oscillator").align_(\center)
.font_(Font("Helvetica", 16));
window.view.decorator.nextLine;

knob_freq = EZKnob(window, 120 @ 90, "freq", ControlSpec(20,440, \lin, step: 10, default: 220),
	{|view| if (sound.notNil) { sound.set(\freq, view.value) }},
	labelWidth: 50);


window.view.decorator.nextLine;


knob_loHarmonics1 = EZKnob(window, 120 @ 90, "loHarmonics1", ControlSpec(0.0, 50.0, \lin, step: 0.25, default: 1.0),
	{|view| if (sound.notNil) { sound.set(\loHarmonics1, view.value) }},
	labelWidth: 50);

knob_hiHarmonics1 = EZKnob(window, 120 @ 90, "hiHarmonics1", ControlSpec(0.0, 50.0, \lin, step: 0.25, default: 6.0),
	{|view| if (sound.notNil) { sound.set(\hiHarmonics1, view.value) }},
	labelWidth: 50);

knob_slope1 = EZKnob(window, 120 @ 90, "slope1", ControlSpec(0.1, 2.0, \lin, step: 0.1, default: 1.0),
	{|view| if (sound.notNil) { sound.set(\slope1, view.value) }},
	labelWidth: 50);

knob_evenOddRatio1 = EZKnob(window, 120 @ 90, "evenOddRatio1", ControlSpec(0.0, 1.0, \lin, step: 0.1, default: 1.0),
	{|view| if (sound.notNil) { sound.set(\evenOddRatio1, view.value) }},
	labelWidth: 50);

knob_spread1 = EZKnob(window, 120 @ 90, "spread1", ControlSpec(1, 5, \lin, step: 1, default: 1),
	{|view| if (sound.notNil) { sound.set(\spread1, view.value) }},
	labelWidth: 50);


window.view.decorator.nextLine;


knob_spreadCompensation = EZKnob(window, 120 @ 90, "spreadCompensation", ControlSpec(0, 1, \lin, step: 1, default: 1),
	{|view| if (sound.notNil) { sound.set(\spreadCompensation, view.value) }},
	labelWidth: 50);

knob_amp = EZKnob(window, 120 @ 90, "amp", ControlSpec(0.1, 1.0, \lin, default: 1.0),
	{|view| if (sound.notNil) { sound.set(\amp, view.value) }},
	labelWidth: 50);

window.view.decorator.left_(10).top_(350);
Button(window, 50@50)
.states_([["Start"], ["Stop", Color.black, Color.red(alpha: 0.2)]])
.action_({|view|  if (view.value ==1 && sound.isNil,
	{sound =
		Synth(\BLOsc1sectionExample,
			[

				\freq, knob_freq.value,

				\loHarmonics1, knob_loHarmonics1.value,
				\hiHarmonics1, knob_hiHarmonics1.value,
				\slope1, knob_slope1.value,
				\evenOddRatio1, knob_evenOddRatio1.value,
				\spread1, knob_spread1.value,

				\spreadCompensation, knob_spreadCompensation.value,

				\mul, knob_amp.value

	])},
	{sound.free; sound = nil }
	)
});
window.front;
)
