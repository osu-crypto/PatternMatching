#pragma once
// This file and the associated implementation has been placed in the public domain, waiving all copyright. No restrictions are placed on its use. 
#include "libOTe/TwoChooseOne/OTExtInterface.h"
#include <cryptoTools/Network/Channel.h>
#include <cryptoTools/Crypto/PRNG.h>
#include <array>

namespace osuCrypto
{

    class PM_IknpOtExtReceiver :
        public OtExtReceiver
    {
    public:
        PM_IknpOtExtReceiver()
            :mHasBase(false)
        {}

        bool hasBaseOts() const override
        {
            return mHasBase;
        }

        bool mHasBase;
        std::array<std::array<PRNG, 2>, gOtExtBaseOtCount> mGens;

        void setBaseOts(
            span<std::array<block, 2>> baseSendOts)override;
        std::unique_ptr<OtExtReceiver> split() override;


        void receive(
            const BitVector& choices,
            span<block> messages,
            PRNG& prng,
            Channel& chl/*,
            std::atomic<u64>& doneIdx*/)override;

    };

}
