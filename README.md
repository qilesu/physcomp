# physcomp
We are trying to tackle the [2018 University Physics Competition](http://www.uphysicsc.com/2018contest.html) Problem B with a python model of a 2D compost pile. We use finite difference Laplacian and Euler's Method to evolve the pile in time, with varying boundary conditions throughout the 24 hour cycle. 

The orginal problem is cited below: 

## Problem B: Compost Pile Sizes
> Composting is the process in which microbes turn organic matter into a useful soil conditioner.  Efficient composting requires the right range of temperatures (40 – 60 degrees C), sufficient moisture content (usually 50-60%), and aeration to deliver adequate oxygen to the microbes.  Larger piles can better maintain high internal temperatures because the heat will diffuse out more slowly.  However, larger piles are also more likely to compact the material, thus inhibiting the flow of oxygen.  What pile size will result in the most efficient composting if we are working with kitchen vegetable waste in a climate where the ambient temperature ranges between 5 and 20 degrees C in a 24-hour period?   How would the most effective pile size vary depending on the ambient temperatures and the input organic materials?

## Our Submission
Our submission can be found in the [submission](/submission) folder. You are welcome to read the paper, and, please don't skip the title page ~~like us~~.

A few highlights are listed here: 

![compressed graph](/submission/compressed-t.png)

The heap actually shrinks down as the bacteria consume it!

![temperature](/submission/long-k.png)

Heat builds up first and dissipates away as the pile shrinks down and bacteria die. 

![bacteria population](/submission/long-o.png)

How the bacteria live and die. 
