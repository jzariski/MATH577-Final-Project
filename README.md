# MATH577-Final-Project

Some important files to note:

coriolanus_encoder1.m: This encoder uses the constructed model from coriolan.txt to arithmetically encode the play as a sequence of states. This method does use hmmtrain but rather recording transitions between states and using that in the econding. Note that in this code we ensure that our trainsition matrix contains no zeros and create beginning estimations for emissions with random numbers. Following the encoding, we test it using a randomly generated sequence using the states from the encoding, and end up with a lossless encoder/decoder that final entropy only slightly above our lower bound of theoretical entropy (different on the order of 1e-2).

corionalus_encoder2.m: This econder is similar to the first one, but has the added step of including hmmtrain with the estimates of transitions and emissions generated in the method from the steps in the first corionalus encoder. Then, using these estimates and hmmtrain, it generates a new estimate for transition probabiliteis and emissions. Then, it encodes a similarlly generated random sequence and losslessly decodes it. Again, the final entropy comes out as only slightly above the lower bound of theoretical entropy ().


