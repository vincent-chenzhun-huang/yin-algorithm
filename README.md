# An Exploration on the YIN Algorithm

See the formatted report [here](https://www.notion.so/An-exploration-on-the-YIN-Algorithm-374119a5d75546df95907f53de1e29ca)

# Installation
- Clone the repository
- Run `pip install -r requirements.txt` to install the dependencies
- Run the Python file yin.py (the whole file or interactively).

# Introduction

With the advancement of different applications that focus extensively on digital signal processing, the ability to react upon the incoming elements has become essential to ensure quality and responsiveness. In music technology industry in particular, frequency / pitch analysis is one of the most important aspects mainly due to the musical nature of various tools such as vocal harmonizer, digital instrument tuner, etc. Unfortunately, the process of developing such tools had been bottlenecked by the accurate detection of the pitch. 

By now, there are plenty of pitch detection libraries available for use in commercial software such as Max/MSP, Matlab. However, there are only a limited number of high-quality implementations out for general public, which makes it non-trivial to build an audio application as simple as a tuner, where all it requires is just an implemented algorithm and some basic UI work. Therefore, I decided to attempt the implementation of one of the most important building blocks of pitch detection, YIN algorithm and gain some insights into different pitch detection algorithms. And in the mean time, test the feasibility of real-time implementation of YIN algorithm.

This report will cover the design of my implementation, optimization techniques, result evaluation, limitations, and conclusions.

# YIN Algorithm

The auto-correlation-based algorithm, YIN was created by Alain de Alain de Cheveigné of IRCAM-CNRS and Hideki Kawahara of Wakayama University. It supports real-time fundamental frequency detection and has a very low error rate as well as very few number of parameters to tune. The creation of this algorithm greatly improves the accuracy of pitch detection and lays solid theorectical foundations for industrial applications that replies extensively on pitch detection algorithms. 

# Project Design

## Technical Stack

Among all the languages that can be used for the implementation, Python and Matlab stood out due to their nature of simple and concise syntax, higher-level control, and fast built-in optimization supports. Those qualities made them suitable for the first time development of the implementation of an algorithm. I chose Python eventually because I was more familiar with it and it gives me more control in both sequential and parallel executions, as opposed to MatLab, which mainly supports vectorized computations. The development environment is Windows Subsystem for Linux (WSL) with Python 3.8.5 64-bit.

## Implementation Steps

The implementation of the YIN algorithm is based on the paper "YIN, a fundamental frequency estimator for speech and music" by Alain de Cheveigne and Hideki Kawahara, and every step of the algorithm replicates the instructions in the paper. There are six steps in total but only the first five steps were implemented to achieve optimal speed in trade of a tiny drop on the precision. The implementation of Step 1 (the auto-correlation method) is omitted in the implementation because Step 2 served as an improvement over this step and does not use the first step. Additionally, during the implementation phase there are underlying preconditions that need to be satisfied for the program to work. `assert` statements were used to ensure the validity of the input arguments. The concept of integration window was also adopted in the paper, so to achieve maximal reuse of the code, the signal is divided into blocks with size of the integration window, and processed by the functions below individually.  Then we combined the processed blocks and form an array of detected frequencies.

The first two steps models the input signal like a period signal and "try out different periods" to find the actual fundamental frequency. We rely on the fact that $x_{t} - x_{t+T} = 0$ for periodic signals and transform this equation into the square sum form

$$ ⁍. $$

Where $\tau$ is the amount of time the signal gets shifted (lag parameter), and W is the size of the integration window. Obviously, this function is equal to 0 at multiples of period and non-zeros elsewhere. We exploit this property and are able to approximately locate the points where we perceive periodicity in an seemingly aperiodic signal. To get the approximated frequency, we simply take the minimum $\tau$ where we get the local minimum of the difference function.  

In the function `diff`, we process one block of signal by comparing it to the next block. Therefore, we cannot have the maximum lag larger than the size of the block processed. Besides, we cannot reach the end of the input signal while processing. So we need 

$$\tau_{max} - 1 + starting\ index + maximum\ index\ in\ current\ window < length\ of\ the\ signal$$

Hence 

$$\tau_{max} < len(x) + 1 - s- w$$

where len(x) is the length of the signal, s is the starting index of the whole input, and w is the integration window size. 

There are scenarios using "raw" difference functions when we erroneously select the zero-lag dip instead of the period dip in the difference function plot and produce incorrect result. The solution proposed in the paper was to replace the difference function with the "cumulative mean normalized difference function" so that the curve starts from 1 and falls below 1 when $d(\tau)$ is below average. This approach eliminates the errors on choosing the starting point of the difference equation and provides a more balanced, and a better-scaled result. The plots of different equations on a single block of signal can be obtained by enabling the `plot` flag on the 

`diff` function. From the example plots, we can see that the "raw" difference function starts from 0 with all the signal points unnormalized while the modified CMNDF starts from 1 as the values with all the signal points scale to 0 - 2.

![An%20exploration%20on%20the%20YIN%20Algorithm%204a2387940a20424583c6ab80144491d0/difference-function-flute-C.png](An%20exploration%20on%20the%20YIN%20Algorithm%204a2387940a20424583c6ab80144491d0/difference-function-flute-C.png)

Difference Function on the last block of a sine wave sample.

![An%20exploration%20on%20the%20YIN%20Algorithm%204a2387940a20424583c6ab80144491d0/cmndf-flute-c.png](An%20exploration%20on%20the%20YIN%20Algorithm%204a2387940a20424583c6ab80144491d0/cmndf-flute-c.png)

CMNDF on the last block of a sine wave sample

It also happens that a higher-order dip of the difference function is deeper than the lower-order dip, which should actually be the detected period. To resolve this issue (octave error), an absolute threshold was introduced in the paper, where we take "the smallest $\tau$ that gives a minimum of $d'$ deeper than that threshold.". This step is equivalent to setting an absolute threshold and take the first $\tau$ that makes the modified difference function fall below that threshold. Therefore, in the implementation, the threshold and `d'` are combined into one variable. 

As the last step of the processing of the block, we adopted to parabolic interpolation to further improve the precision of the pitch detection algorithm. In each block of signal, we fit the local minimum and its immediate neighbors in a parabola. Then we used the interpolated ordinates to get the final result. Note that the last step has been omitted in the implementation for the sake of efficiency. And the trade-off is only a slight drop in accuracy.

Finally, the processing results from each block is combined to form an array of detected frequencies. Since WSL does not support audio I/O, I could not produce sine waves with specified frequencies in place of the fundamental frequencies detected. Instead, a plot with the detected frequencies vs block is created to perform a preliminary evaluation on the results of the algorithm.

# Results and Evaluations

### Pitch Detection on Different Samples

For the testing sound samples, I obtained clean, processed samples from [http://wiki.laptop.org/go/Sound_samples](http://wiki.laptop.org/go/Sound_samples) and use different files as the input to the implemented algorithm as a preliminary evaluation. Here are some selected plots coming out of the script.

![An%20exploration%20on%20the%20YIN%20Algorithm%204a2387940a20424583c6ab80144491d0/flute-C.png](An%20exploration%20on%20the%20YIN%20Algorithm%204a2387940a20424583c6ab80144491d0/flute-C.png)

Pitch detected with flute playing C

![An%20exploration%20on%20the%20YIN%20Algorithm%204a2387940a20424583c6ab80144491d0/sin-wave.png](An%20exploration%20on%20the%20YIN%20Algorithm%204a2387940a20424583c6ab80144491d0/sin-wave.png)

Pitch detected with the sum of three sine waves at 400Hz, 800Hz and 1200Hz

![An%20exploration%20on%20the%20YIN%20Algorithm%204a2387940a20424583c6ab80144491d0/shenai-C.png](An%20exploration%20on%20the%20YIN%20Algorithm%204a2387940a20424583c6ab80144491d0/shenai-C.png)

Pitch detected with shenai playing C

![An%20exploration%20on%20the%20YIN%20Algorithm%204a2387940a20424583c6ab80144491d0/trumpet-lick.png](An%20exploration%20on%20the%20YIN%20Algorithm%204a2387940a20424583c6ab80144491d0/trumpet-lick.png)

Pitches detected with a trumpet lick

As we can see from the graph, the detected pitches for the first three WAV files are very close to the actual fundamental frequencies generated within (10%). And the output with trumpet lick shows a decent range of detected pitches with one “incident" around chunk 2, where the frequency just dropped drastically. However, overall, the plot sketched out the progression of the lick and was able to detect most (if not all) notes in the audio file. For the generated sine wave, in particular, the detected frequency stabilized at 405 Hz, which implies the validity of the algorithm on periodic signal inputs.

However, the implementation does not work for unprocessed, noisy samples. For example, it was not able to detect the pitch from a piano playing C note in a noisy background. Moreover, this implementation does not work on the host environment (Windows 10) for some reason. When I ran it on windows with the exact same configurations, I got noisy output. I have not had the chance to test it on the Mac systems yet but it should work the same way as it does on WSL.

### Runtime and Optimizations

Since YIN is an algorithm that is designed for real-time performance, testing its runtime is necessary to make sure that the delay is acceptable from a theoretical point of view. 

In the implementation, I use an audio block of 4410 samples to test the execution time of the algorithm. Assume that we are working with real-time audio I/O streams and the sample rate is 44100Hz, the amount of time taken to obtain 4410 samples is

$$ ⁍. $$

With Python's built-in module `time`, we can easily get the execution time of the implemented YIN algorithm. If everything is executed sequentially, the runtime is 0.4s testing on flute sample for every 0.1s fragment. This is, of course, unacceptable. 

Note that Python supports vectorizations with Numpy, where we could parallelize computations that do not affect each other. There is a big bottleneck on the runtime while calculating the difference function between the current and the next block because everything is sequential when they do not affect each other. The solution is to put them in NumPy arrays and parallelize the application of the difference function. A similar approach is used on `cmndiff` . Vectorization improved the runtime to 0.06s per 0.1s sample. This result gives the implementation great potential for real-time performance if we can parallelize I/O and computations with some care.

If real-time performance is not required, further optimization can also be achieved with the utilization of multiple threads. The current implementation only uses one main thread to divide, compute and combine all the blocks in the signal, which can be a drag on the performance if we have a very long signal input. A potential solution is to divide the read input signal and spawn a number of threads to apply parallel block-wise processing, then combine the predicted results into one final array. This approach can be faster than using only vectorization but it creates possibility for data race on the shared array. Locking could potentially even slows down the execution if we have a lot of data fragments. Besides, Numpy already supports parallel computing with threads and it could cause problems when working with a black box like that. For all the risks above and limit on time, I decided not to implement this optimization.

# Conclusion

Overall, this implementation could be used as a testimony that the YIN algorithm does work well on certain inputs. This attempt has also proved the implementation of pitch detection algorithms non-trivial. The implementation can be platform-dependent and we have to interact with hardware devices to create audio I/O. Moreover, a lot of care on the optimization techniques and analysis of the variable validity is required for the algorithm to work properly. 

Future work on this implementation can focus on cross-platform support and reusability. Famous frameworks like JUCE does not provide pitch tracking implementations, so this can also be a potential addition besides adding third-party libraries like STK.
