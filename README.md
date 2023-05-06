# MATH577-Final-Project

Some important files to note:

coriolanus_encoder1.m: This encoder uses the constructed model from coriolan.txt to arithmetically encode the play as a sequence of states. This method does not use hmmtrain but rather records transitions between states and uses that in the econding. Note that in this code we ensure that our trainsition matrix contains no zeros and create beginning estimations for emissions with random numbers. Following the encoding, we test it using a randomly generated sequence using the states from the encoding, and end up with a lossless encoder/decoder with final entropy only slightly above our lower bound of theoretical entropy (different on the order of 1e-2).

corionalus_encoder2.m: This econder is similar to the first one, but has the added step of including hmmtrain with the estimates of transitions and emissions generated in the method from the steps in the first corionalus encoder. Then, using these estimates and hmmtrain, it generates a new estimate for transition probabiliteis and emissions. Then, it encodes a similarlly generated random sequence and losslessly decodes it. Again, the final entropy comes out as only slightly above the lower bound of theoretical entropy (also different on the order of 1e-2). However, the difference was still less than half that of the first method, so I can conclude that this method is relatively more efficient. 


After testing both of these methods, I found that both methods had an entropy difference on the same order of magnitude, though the second encoder/decoder method was slightly more efficient. Both were very close to the theoretical entropy value, much more so than regular arithmetic coding, so I would say that both encoding methods are fairly convenient. Note that in both cases I stuck with k=2, since it allowed for the code to run faster then higher sizes of states. I also limited the iterations of hmmtrain, and I would expect increasing said iterations would add efficiency to the second encoder.

In this github I included the provided cicero and script1 methods detailing the encoding of the HMM. In the coriolanus_encoder1.m and coriolanus_encoder2.m files, I commented in places where I made significant changes, such as altering the count number in decoding states when it (rarely) defaults to zero. 
