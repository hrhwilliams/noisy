#ifndef __NOISE_MASK_HPP_
#define __NOISE_MASK_HPP_

namespace Noise {
    class Mask {
    protected:
    public:
        Mask();
    };

    class MaskingBox : public Mask {
    private:
    public:
        MaskingBox();
    };
};

#endif // __NOISE_MASK_HPP_
