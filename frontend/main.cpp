
//using namespace std;
#include "main.h"
#include <chrono>
#include <thread>

using namespace std;
#include <fstream>

//#define PRINT
#define TIME
//#define TEST

//#########NOTE: the below values effect on the perfomance. Depend on bandwith, need to manually change it

u64 stepSize = 1 << 12; //LAN: 1<<12; WAN: 1<<20
u64  blk_step = 4; //LAN: 4; WAN: 
u64  stepDeltaSize = 4; //LAN: 4; WAN: 16

u64 mMaskSize = 40 / 8;
u64 mNumBlkText = 1<< (14 - 7);// //1563 (4x100,000/128/2) ,14063 (36x100,000/128/2)
u64 mNumBlkPattern = 1<< (9 - 7); //32
u64 numTrials =1;

void pmSender(BitVector choicesTextOn, int numBlkText, int numBlkPattern)
{
	for (u64 idTrial = 0; idTrial < numTrials; idTrial++)
	{
		std::this_thread::sleep_for(std::chrono::seconds(3));

#pragma region channel
		setThreadName("Sender");
		IOService ios(0);
		Endpoint ep1(ios, "127.0.0.1", 1212, EpMode::Client, "ep");
		Channel senderChannel = ep1.addChannel("chl", "chl");
#pragma endregion

		Timer mTimer;
		PRNG prng0(_mm_set_epi32(4253465, 3434565, 234435, 23987045));
		PRNG prng1(_mm_set_epi32(4253233465, 334565, 0, 235));

		u64 numOTs = numBlkPattern * 128;
		u64 numBitText = numBlkText * 128;
		int numOPPRFs = numBitText - numOTs + 1;
		int numBlkOPPRFs = numBlkText - numBlkPattern + 1;





		std::vector<block> recvSeeds(numOTs), baseOTRecv(128);
		std::vector<std::array<block, 2>> sendSeeds(numOTs), baseOTSend(128);

		block* pmSendMsg = new block[numOTs*numBlkOPPRFs];
		block* trans_pmSendMsg = new block[numBlkOPPRFs*numOTs];
		block* blkTextOff = new block[numBlkText];
		block* blkTextOn = new block[numBlkText];


		BitVector baseOTChoice(128),
			choicesTextOff(numBitText);

		baseOTChoice.randomize(prng0);
		choicesTextOff.randomize(prng0);

#ifdef TEST
		for (u64 i = 0; i < numOTs; i++)
		{
			choicesTextOff[i] = choicesTextOn[i];
		}

#endif // TEST

		for (u64 i = 0; i < numBlkText; i++)
		{
			blkTextOff[i] = toBlock(choicesTextOff.data() + (i * sizeof(block)));
			blkTextOn[i] = toBlock(choicesTextOn.data() + (i * sizeof(block)));
		}


		IknpOtExtSender senderOT;
		KkrtNcoOtSender kkrtSender;

		kkrtSender.configure(false, 40, 128);

		// the number of base OT that need to be done
		u64 baseOPRFCount = kkrtSender.getBaseOTCount();

		std::vector<block> baseOPRFRecv(baseOPRFCount);
		BitVector baseOPRFChoice(baseOPRFCount);
		baseOPRFChoice.randomize(prng0);


		NaorPinkas baseOTs;

		baseOTs.receive(baseOTChoice, baseOTRecv, prng0, senderChannel, 1);
		senderOT.setBaseOts(baseOTRecv, baseOTChoice);
		senderOT.send(sendSeeds, prng1, senderChannel);

		//base-oprf
		std::array<std::array<block, 2>, 128> baseBaseOT;
		senderOT.send(baseBaseOT, prng1, senderChannel);

		IknpOtExtReceiver base;
		base.setBaseOts(baseBaseOT);
		base.receive(baseOPRFChoice, baseOPRFRecv, prng1, senderChannel);

		//oprf
		kkrtSender.setBaseOts(baseOPRFRecv, baseOPRFChoice);
		kkrtSender.init(numOPPRFs, prng1, senderChannel);

		//matrix |p|x|t|
		u64 idRow = 0;
		std::vector<PRNG> prg_pm(2);
		std::vector<block*> genS(2);
		genS[0] = new block[numBlkOPPRFs];
		genS[1] = new block[numBlkOPPRFs];

		vector<block*> textOffShift(128);
		for (u64 j = 0; j < 128; j += 1)
		{
			textOffShift[j] = new block[numBlkText];
			blks_bitshift_left(blkTextOff, textOffShift[j], numBlkText, j); //blkTextOn[1...p]			

		}


		while (idRow < numOTs)
		{
			u64 currStepDeltaSize = std::min(stepDeltaSize, numOTs - idRow);

			//block* delta_send = new block[numBlkText*currStepDeltaSize];

			uPtr<Buff> sendDeltaBuff(new Buff);

			sendDeltaBuff->resize(((numBlkOPPRFs - 1) * sizeof(block) + sizeof(u8))*currStepDeltaSize);
			auto deltaView = sendDeltaBuff->getMatrixView<u8>((numBlkOPPRFs - 1) * sizeof(block) + sizeof(u8));

			for (u64 j = 0; j < currStepDeltaSize; j++)
			{
				prg_pm[0].SetSeed(sendSeeds[idRow][0]);
				prg_pm[1].SetSeed(sendSeeds[idRow][1]);


				prg_pm[0].get(genS[0], numBlkOPPRFs);
				prg_pm[1].get(genS[1], numBlkOPPRFs);

				for (u64 k = 0; k < numBlkOPPRFs-1; ++k)
				{
					pmSendMsg[idRow*numBlkOPPRFs + k] = genS[0][k] ^ textOffShift[idRow % 128][idRow / 128 + k];
					block tmp = pmSendMsg[idRow*numBlkOPPRFs + k] ^ genS[1][k];

					
						memcpy(deltaView[j].data() + k * sizeof(block), &tmp, sizeof(block));

				}

				u64 k = numBlkOPPRFs - 1;
				pmSendMsg[idRow*numBlkOPPRFs + k] = genS[0][k] ^ textOffShift[idRow % 128][idRow / 128 + k];
				block tmp = pmSendMsg[idRow*numBlkOPPRFs + k] ^ genS[1][k];
				memcpy(deltaView[j].data() + k * sizeof(u8), &tmp, sizeof(u8));



				//shift
				//blks_bitshift_left(pmSendMsg + idRow*numBlkText, numBlkText, idRow);

				idRow++;
			}
			//	cout <<numOTs << ": " << idRow << " ss " << currStepDeltaSize << "\n";
			//senderChannel.send((u8*)delta_send, sizeof(block)*(numBlkText*currStepDeltaSize));
			senderChannel.asyncSend(std::move(sendDeltaBuff));
			//cout << numOTs << " ss1 " << currStepDeltaSize << "\n";
		}
		
		//trans
		u64 j = 0;
		while (j < numBlkOPPRFs)
		{
			u64 cur_blk_step = std::min(blk_step, numBlkOPPRFs - j);

			block* tmp_s = new block[numOTs*cur_blk_step];
			for (u64 i = 0; i < numOTs; i++)
				memcpy(tmp_s + i*cur_blk_step, pmSendMsg + i*numBlkOPPRFs + j, cur_blk_step * sizeof(block));

			sse_trans((uint8_t*)(trans_pmSendMsg + j*numOTs), (uint8_t*)tmp_s, numOTs, 128 * cur_blk_step);

			j += cur_blk_step;
		}
		delete[] pmSendMsg;


#ifdef TEST

		senderChannel.send((u8*)trans_pmSendMsg, sizeof(block)*(numBlkOPPRFs*numOTs));


#endif
		//####################
		//#####Online Phase###
		//####################


		block* corrText = new block[numBlkText];//blkTextOn online
		block* corrWild = new block[numBlkPattern];

		for (u64 j = 0; j < numBlkText; j += 1)
			corrText[j] = blkTextOff[j] ^ blkTextOn[j];

		senderChannel.send((u8*)corrText, sizeof(block)*(numBlkText));


#ifdef TEST
		std::cout << IoStream::lock;
		cout << "s corrText " << corrText[numBlkText - 1] << "\n";
		std::cout << IoStream::unlock;

#endif

		//for (u64 j = 0; j < numOPPRFs; j += 1)
		//{
		//	if (j != 0)
		//		blks_bitshift_left(corrText, numBlkText, 1); //blkTextOn[1...p]			


		//	for (u64 i = 0; i < numBlkPattern; i += 1)
		//		trans_pmSendMsg[j*numBlkPattern + i] = trans_pmSendMsg[j*numBlkPattern + i] ^ corrText[i];
		//}

		//vector<block*> corrTextShift(128);
		for (u64 j = 0; j < 128; j += 1)
		{
			//corrTextShift[j] = new block[numBlkText];
			blks_bitshift_left(corrText, textOffShift[j], numBlkText, j); //blkTextOn[1...p]			

		}

		for (u64 j = 0; j < numOPPRFs; j += 1)
		{
			for (u64 i = 0; i < numBlkPattern; i += 1)
			{
				trans_pmSendMsg[j*numBlkPattern + i] = trans_pmSendMsg[j*numBlkPattern + i] ^ textOffShift[j % 128][j / 128 + i];
			}
		}

		senderChannel.recv((u8*)corrWild, sizeof(block)*(numBlkPattern));

#ifdef TEST
		std::cout << IoStream::lock;
		cout << "s corrWild " << corrWild[numBlkPattern - 1] << "\n";
		std::cout << IoStream::unlock;

#endif


		//for (u64 j = 0; j < numOPPRFs; j += 1)
		//{
		//	if (j != 0)
		//		blks_bitshift_left(blkTextOn, numBlkText, 1); //blkTextOn[1...p]	

		//	for (u64 i = 0; i < numBlkPattern; i += 1)
		//		trans_pmSendMsg[j*numBlkPattern + i] = trans_pmSendMsg
		//		[j*numBlkPattern + i] ^ (corrWild[i] & blkTextOn[i]);
		//}

		//vector<block*> textShift(128);
		for (u64 j = 0; j < 128; j += 1)
		{
		//	textShift[j] = new block[numBlkText];
			blks_bitshift_left(blkTextOn, textOffShift[j], numBlkText, j); //blkTextOn[1...p]			

		}

		for (u64 j = 0; j < numOPPRFs; j += 1)
		{
			for (u64 i = 0; i < numBlkPattern; i += 1)
			{
				trans_pmSendMsg[j*numBlkPattern + i] = trans_pmSendMsg
					[j*numBlkPattern + i] ^ (corrWild[i] & textOffShift[j % 128][j / 128 + i]);
			}
		}




#pragma region OPRF

		//####################
		//#####OPRF###
		//####################

#if 1



#if 1
		static const int oprfWidth(4);

		u64 i = 0;
		while (i<numOPPRFs)
			//for (u64 i = 0; i < numOPPRFs; i += stepSize)
		{
			u64 currStepSize = std::min(stepSize, numOPPRFs - i);

			block* encoding2 = new block[currStepSize];
			std::vector<block*> sendInputs(currStepSize);
			uPtr<Buff> sendMaskBuff(new Buff);

			sendMaskBuff->resize(currStepSize  * mMaskSize * sizeof(u8));
			auto maskView = sendMaskBuff->getMatrixView<u8>(mMaskSize);



			for (u64 j = 0; j < currStepSize; j++)
			{
				sendInputs[j] = new block[oprfWidth];

				if (numBlkPattern >= oprfWidth)
				{ //copy first (oprfWidth-1) pmMsg to OPRF
					memcpy((u8*)sendInputs[j]
						, (u8*)(trans_pmSendMsg + (i + j)*numBlkPattern), (oprfWidth) * sizeof(block));

					//compute XOR the rest
					for (u64 k = oprfWidth; k < numBlkPattern; k++)
					{
						sendInputs[j][oprfWidth - 1] = sendInputs[j][oprfWidth - 1]
							^ trans_pmSendMsg[(i + j)*numBlkPattern + k];

					}
				}
				else
				{
					//copy first pmMsg to OPRF
					memcpy((u8*)sendInputs[j]
						, (u8*)(trans_pmSendMsg + (i + j)*numBlkPattern), (numBlkPattern) * sizeof(block));

					//compute AES for the rest
					block* tmp_send = new block[oprfWidth - numBlkPattern];

					block tmp_send_plain = ZeroBlock;
					for (u64 k = 0; k < numBlkPattern; k++)
					{
						tmp_send_plain = tmp_send_plain ^ sendInputs[j][k];

					}

					std::array<block, oprfWidth> choice_send{ tmp_send_plain,tmp_send_plain ,tmp_send_plain ,tmp_send_plain };
					kkrtSender.mMultiKeyAES.ecbEncNBlocks(choice_send.data(), sendInputs[j]);
				}
			}


			// receive the next stepSize correction values. This allows the sender to now call encode
			// on the next stepSize OTs.
			kkrtSender.recvCorrection(senderChannel, currStepSize);


			for (u64 k = 0; k < currStepSize; ++k)
			{
				// the sender can now call encode(i, ...) for k \in {0, ..., i}. 
				// Lets encode the same input and then we should expect to
				// get the same encoding.

				kkrtSender.pmEncode(i + k, sendInputs[k], maskView[k].data(), mMaskSize);
			}

			senderChannel.asyncSend(std::move(sendMaskBuff));
			i += currStepSize;
		}


#endif



#endif
#pragma endregion

#ifdef TEST
		//for test
		senderChannel.send((u8*)trans_pmSendMsg, sizeof(block)*(numBlkOPPRFs*numOTs));
#endif


		delete[] trans_pmSendMsg;
		delete[] blkTextOff;
		delete[] blkTextOn;
		for (u64 j = 0; j < 128; j += 1)
		{
			delete[] textOffShift[j];
		}

		senderChannel.close(); ep1.stop();  ios.stop();
	}
}

void pmRecv(BitVector choicesPatternOn, BitVector choicesWildOn, int numBlkText, int numBlkPattern)
{
	double avg_total_time_ms = 0;
	double avg_ontime_ms = 0;
	std::fstream fOutputs, fLatex;
	fOutputs.open("./output.txt", fOutputs.app | fOutputs.out);
	fLatex.open("./fLatex.txt", fLatex.app | fLatex.out);

	for (u64 idTrial = 0; idTrial < numTrials; idTrial++)
	{
		std::this_thread::sleep_for(std::chrono::seconds(3));


#pragma region channel
		setThreadName("Sender");
		IOService ios(0);
		Endpoint ep0(ios, "127.0.0.1", 1212, EpMode::Server, "ep");
		Channel recvChannel = ep0.addChannel("chl", "chl");
#pragma endregion

		Timer mTimer;

		PRNG prng0(_mm_set_epi32(4253465, 3434565, 234435, 23987045));
		PRNG prng1(_mm_set_epi32(4253233465, 334565, 0, 235));

		u64 numOTs = numBlkPattern * 128;
		u64 numBitText = numBlkText * 128;
		int numOPPRFs = numBitText - numOTs + 1;

		int numBlkOPPRFs = numBlkText - numBlkPattern + 1;



		std::vector<block> recvSeeds(numOTs), baseOTRecv(128);
		std::vector<std::array<block, 2>> sendSeeds(numOTs), baseOTSend(128);

		block* pmRecvMsg = new block[numOTs*numBlkOPPRFs];
		block* trans_pmRecvMsg = new block[numBlkOPPRFs*numOTs];
		block* trans_pmSendMsg = new block[numBlkOPPRFs*numOTs];
		block* blkPatternOff = new block[numBlkPattern];
		block* blkPatternOn = new block[numBlkPattern];
		block* blkWildOff = new block[numBlkPattern];
		block* blkWildOn = new block[numBlkPattern];


		BitVector baseOTChoice(128),
			choicesWildOff(numOTs),
			choicesPatternOff(numOTs);

		baseOTChoice.randomize(prng0);
		choicesWildOff.randomize(prng0);
		choicesPatternOff.randomize(prng0);

#ifdef TEST
		for (u64 i = 0; i < numOTs; i++)
		{
			choicesWildOff[i] = choicesWildOn[i];
			choicesPatternOff[i] = choicesPatternOn[i];
		}

#endif // TEST



		for (u64 i = 0; i < numBlkPattern; i++)
		{
			blkPatternOff[i] = toBlock(choicesPatternOff.data() + (i * sizeof(block)));
			blkWildOff[i] = toBlock(choicesWildOff.data() + (i * sizeof(block)));
			blkPatternOn[i] = toBlock(choicesPatternOn.data() + (i * sizeof(block)));
			blkWildOn[i] = toBlock(choicesWildOn.data() + (i * sizeof(block)));
		}



		IknpOtExtReceiver recvOT;
		KkrtNcoOtReceiver kkrtRecv;

		kkrtRecv.configure(false, 40, 128);


		// the number of base OT that need to be done
		u64 baseOPRFCount = kkrtRecv.getBaseOTCount();

		// Fake some base OTs
		std::vector<block> baseOPRFRecv(baseOPRFCount);
		std::vector<std::array<block, 2>> baseOPRFSend(baseOPRFCount);
		BitVector baseOPRFChoice(baseOPRFCount);
		baseOPRFChoice.randomize(prng0);

		NaorPinkas baseOTs;

#ifdef TIME 
		auto start = mTimer.setTimePoint("offline");
#endif

		baseOTs.send(baseOTSend, prng1, recvChannel, 1);
		recvOT.setBaseOts(baseOTSend);
		recvOT.receive(choicesWildOff, recvSeeds, prng0, recvChannel);

		////base oprf
		std::array<block, 128> baseBaseOT;
		BitVector baseBaseChoice(128);
		baseBaseChoice.randomize(prng1);
		recvOT.receive(baseBaseChoice, baseBaseOT, prng1, recvChannel);

		IknpOtExtSender base;
		base.setBaseOts(baseBaseOT, baseBaseChoice);
		base.send(baseOPRFSend, prng1, recvChannel);

		//oprf
		kkrtRecv.setBaseOts(baseOPRFSend);
		kkrtRecv.init(numOPPRFs, prng1, recvChannel);


		//matrix |p|x|t|
#ifdef TIME
		auto startShift = mTimer.setTimePoint("startShift");
#endif
		u64 idRow = 0;
		ByteStream deltaBuffer;
		PRNG prg_pm;

		while (idRow < numOTs)
		{
			u64 currStepDeltaSize = std::min(stepDeltaSize, numOTs - idRow);

			recvChannel.recv(deltaBuffer);
			auto deltaView = deltaBuffer.getMatrixView<u8>((numBlkOPPRFs - 1) * sizeof(block) + sizeof(u8));


			for (u64 j = 0; j < currStepDeltaSize; j++)
			{
				prg_pm.SetSeed(recvSeeds[idRow]);
				prg_pm.get(pmRecvMsg + idRow*numBlkOPPRFs, numBlkOPPRFs);

				block tmp = ZeroBlock;

				if (choicesWildOff[idRow] == 1)
				{
					for (u64 k = 0; k < numBlkOPPRFs-1; ++k)
					{
						memcpy(&tmp, deltaView[j].data() + k * sizeof(block), sizeof(block));
						pmRecvMsg[idRow*numBlkOPPRFs + k] = pmRecvMsg[idRow*numBlkOPPRFs + k] ^ tmp;
					}

					u64 k = numBlkOPPRFs - 1;
					memcpy(&tmp, deltaView[j].data() + k * sizeof(u8), sizeof(u8));
					pmRecvMsg[idRow*numBlkOPPRFs + k] = pmRecvMsg[idRow*numBlkOPPRFs + k] ^ tmp;

				}
				else
					if (choicesPatternOff[idRow] == 1)
						for (u64 k = 0; k < numBlkOPPRFs; ++k)
							pmRecvMsg[idRow*numBlkOPPRFs + k] = pmRecvMsg[idRow*numBlkOPPRFs + k] ^ AllOneBlock;

				//blks_bitshift_left(pmRecvMsg + idRow*numBlkText, numBlkText, idRow);

				idRow++;
			}
		}

#ifdef TIME 
		auto endShift = mTimer.setTimePoint("endShift");
		auto timeTmp = std::chrono::duration_cast<std::chrono::milliseconds>(endShift - startShift).count();
		std::cout << "endShift= " << timeTmp << "ms\n";
#endif

		auto startTrans = mTimer.setTimePoint("startTrans");
		//trans
		u64 j = 0;
		while (j < numBlkOPPRFs)
		{
			u64 cur_blk_step = std::min(blk_step, numBlkOPPRFs - j);

			block* tmp_s = new block[numOTs*cur_blk_step];

			for (u64 i = 0; i < numOTs; i++)
				memcpy(tmp_s + i*cur_blk_step, pmRecvMsg + i*numBlkOPPRFs + j, cur_blk_step * sizeof(block));

			sse_trans((uint8_t*)(trans_pmRecvMsg + j*numOTs), (uint8_t*)tmp_s, numOTs, 128 * cur_blk_step);

			j += cur_blk_step;
		}

		delete[] pmRecvMsg;

#ifdef TIME 
		auto endTrans = mTimer.setTimePoint("endTrans");
		timeTmp = std::chrono::duration_cast<std::chrono::milliseconds>(endTrans - startTrans).count();
		std::cout << "endTrans= " << timeTmp << "ms\n";
#endif

#ifdef TEST
		//block* trans_pmSendMsg = new block[numOTs*numBlkText];

		recvChannel.recv((u8*)trans_pmSendMsg, sizeof(block)*(numBlkOPPRFs*numOTs));


		for (u64 i = 0; i < numOPPRFs; ++i)
			for (u64 j = 0; j < numBlkPattern; ++j)
			{
				block revcBlock = trans_pmRecvMsg[i*numBlkPattern + j];
				block senderBlock = trans_pmSendMsg[i*numBlkPattern + j];

				if (eq(revcBlock, senderBlock))
					cout << "off_trans-- " << i << " " << j << ": " << senderBlock << "  " << revcBlock << endl;

			}

		u64 ii = 0; u64 jj = 7;
		//if (eq(revcBlock, senderBlock))
		cout << "off_trans-- " << ii << " " << jj << ": " << trans_pmSendMsg[ii*numBlkPattern + jj] << "  "
			<< trans_pmRecvMsg[ii*numBlkPattern + jj] << endl;

#endif

		//####################
		//#####Online Phase###
		//####################

#ifdef TIME 
		auto startInputOn = mTimer.setTimePoint("startInputOn");
#endif

		block* corrText = new block[numBlkText];
		block* corrPattern = new block[numBlkPattern];
		block* corrWild = new block[numBlkPattern];

		for (u64 j = 0; j < numBlkPattern; j += 1)
			corrPattern[j] = blkPatternOn[j] ^ blkPatternOff[j];

		for (u64 j = 0; j < numBlkPattern; j += 1)
			corrWild[j] = blkWildOn[j] ^ blkWildOff[j];

		recvChannel.recv((u8*)corrText, sizeof(block)*(numBlkText));

#ifdef TEST
		std::cout << IoStream::lock;
		cout << "r corrText " << corrText[numBlkText - 1] << "\n";
		cout << "r corrWild " << corrWild[numBlkPattern - 1] << "\n";
		std::cout << IoStream::unlock;
#endif

#ifdef TIME 
		auto startTextOn = mTimer.setTimePoint("startTextOn");
		timeTmp = std::chrono::duration_cast<std::chrono::milliseconds>(startTextOn - startInputOn).count();
		std::cout << "startTextOn - startInputOn= " << timeTmp << "ms\n";
#endif

		vector<block*> textShift(128);
		for (u64 j = 0; j < 128; j += 1)
		{
			textShift[j] = new block[numBlkText];
			blks_bitshift_left(corrText, textShift[j], numBlkText, j); //blkTextOn[1...p]			

		}

		for (u64 j = 0; j < numOPPRFs; j += 1)
		{
			for (u64 i = 0; i < numBlkPattern; i += 1)
			{
				trans_pmRecvMsg[j*numBlkPattern + i] = trans_pmRecvMsg
					[j*numBlkPattern + i] ^ (blkWildOff[i] & textShift[j % 128][j / 128 + i]);
			}
		}

#ifdef TIME 
		auto endTextOn = mTimer.setTimePoint("endTextOn");
		timeTmp = std::chrono::duration_cast<std::chrono::milliseconds>(endTextOn - startTextOn).count();
		std::cout << "endTextOn= " << timeTmp << "ms\n";
#endif

		//pattern on: 
		//receiver locally flips the bit where patternOff!=patternOn and !wild
		for (u64 j = 0; j < numOPPRFs; j += 1)
			for (u64 i = 0; i < numBlkPattern; i += 1)
				trans_pmRecvMsg[j*numBlkPattern + i] = trans_pmRecvMsg
				[j*numBlkPattern + i] ^ ((blkWildOff[i] ^ AllOneBlock)&corrPattern[i]);

#ifdef TIME 
		auto endPattternOn = mTimer.setTimePoint("endPattternOn");
		timeTmp = std::chrono::duration_cast<std::chrono::milliseconds>(endPattternOn - endTextOn).count();
		std::cout << "endPattternOn= " << timeTmp << "ms\n";
#endif


		recvChannel.send((u8*)corrWild, sizeof(block)*(numBlkPattern));

		for (u64 j = 0; j < numOPPRFs; j += 1)
			for (u64 i = 0; i < numBlkPattern; i += 1)
				trans_pmRecvMsg[j*numBlkPattern + i] = trans_pmRecvMsg
				[j*numBlkPattern + i] ^ ((corrWild[i] & blkPatternOn[i]));

#ifdef TIME 
		auto endWilOn = mTimer.setTimePoint("endWilOn");
		timeTmp = std::chrono::duration_cast<std::chrono::milliseconds>(endWilOn - endPattternOn).count();
		std::cout << "endWilOn= " << timeTmp << "ms\n";

		auto endInputOn = mTimer.setTimePoint("endInputOn");
		timeTmp = std::chrono::duration_cast<std::chrono::milliseconds>(endInputOn - startInputOn).count();
		std::cout << "endInputOn= " << timeTmp << "ms\n";

#endif




#pragma region OPRF

		//####################
		//#####OPRF###
		//####################

#ifdef TIME 
		auto startOprf = mTimer.setTimePoint("endOprf");
#endif
		static const int oprfWidth(4);
		cout << "numOPPRFs: " << numOPPRFs << endl;

		MatrixView<u8> maskView;
		ByteStream maskBuffer;

		u64 i = 0;

		while (i < numOPPRFs)
			//	for (u64 i = 0; i < numOPPRFs; i += stepSize)
		{
			auto startOprfInput = mTimer.setTimePoint("startOprfInput");
			u64 currStepSize = std::min(stepSize, numOPPRFs - i);
			std::vector<block*> recvInputs(currStepSize);
			block* encoding1 = new block[currStepSize];


			for (u64 j = 0; j < currStepSize; j++)
			{
				recvInputs[j] = new block[oprfWidth];

				if (numBlkPattern >= oprfWidth)
				{ //copy first (oprfWidth-1) pmMsg to OPRF

					memcpy((u8*)recvInputs[j]
						, (u8*)(trans_pmRecvMsg + (i + j)*numBlkPattern), (oprfWidth) * sizeof(block));

					//compute XOR the rest
					for (u64 k = oprfWidth; k < numBlkPattern; k++)
					{
						recvInputs[j][oprfWidth - 1] = recvInputs[j][oprfWidth - 1]
							^ trans_pmRecvMsg[(i + j)*numBlkPattern + k];
					}
				}
				else
				{
					//copy first pmMsg to OPRF
					memcpy((u8*)recvInputs[j]
						, (u8*)(trans_pmRecvMsg + (i + j)*numBlkPattern), (numBlkPattern) * sizeof(block));

					//compute AES for the rest
					block* tmp_recv = new block[oprfWidth - numBlkPattern];
					block tmp_recv_plain = ZeroBlock;
					for (u64 k = 0; k < numBlkPattern; k++)
					{
						tmp_recv_plain = tmp_recv_plain ^ recvInputs[j][k];
					}

					std::array<block, oprfWidth> choice_recv{ tmp_recv_plain,tmp_recv_plain ,tmp_recv_plain ,tmp_recv_plain };
					kkrtRecv.mMultiKeyAES.ecbEncNBlocks(choice_recv.data(), recvInputs[j]);
				}
			}

			auto endOprfInput = mTimer.setTimePoint("endOprfInput");
			/*timeTmp = std::chrono::duration_cast<std::chrono::milliseconds>(endOprfInput - startOprfInput).count();
			std::cout << "endOprfInput= " << timeTmp << "ms\n";*/


			for (u64 k = 0; k < currStepSize; ++k)
			{

				// The receiver MUST encode before the sender. Here we are only calling encode(...) 
				// for a single i. But the receiver can also encode many i, but should only make one 
				// call to encode for any given value of i.

				kkrtRecv.pmEncode(i + k, recvInputs[k], (u8*)&encoding1[k], mMaskSize);
			}

			// This call will send to the other party the next "stepSize" corrections to the sender.
			// If we had made more or less calls to encode above (for contigious i), then we should replace
			// stepSize with however many calls we made. In an extreme case, the reciever can perform
			// encode for i \in {0, ..., numOTs - 1}  and then call sendCorrection(recvChl, numOTs).
			kkrtRecv.sendCorrection(recvChannel, currStepSize);



			recvChannel.recv(maskBuffer);
			maskView = maskBuffer.getMatrixView<u8>(mMaskSize);

			//// check that we do in fact get the same value
			for (u64 k = 0; k < currStepSize; ++k)
				if (memcmp(maskView[k].data(), &encoding1[k], mMaskSize) == 0)
					//	if (eq(encoding1r[k], encoding1[k]))
				{
					cout << "ouput matching oprf-- " << i + k << ": " << encoding1[k] << endl;

				}

			i += currStepSize;
		}

#ifdef TIME
		auto endOprf = mTimer.setTimePoint("endOprf");
		timeTmp = std::chrono::duration_cast<std::chrono::milliseconds>(endOprf - startOprf).count();
		std::cout << "endOprf= " << timeTmp << "ms\n";
#endif

#pragma endregion

#ifdef TIME 
		double time1 = std::chrono::duration_cast<std::chrono::milliseconds>(endTrans - start).count();
		double time2 = std::chrono::duration_cast<std::chrono::milliseconds>(endInputOn - startInputOn).count();
		double time3 = std::chrono::duration_cast<std::chrono::milliseconds>(endOprf - startOprf).count();
		double time_second = std::chrono::duration_cast<std::chrono::seconds>(endOprf - start).count();
		double ontime_second = std::chrono::duration_cast<std::chrono::seconds>(endOprf - startInputOn).count();

		double time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(endOprf - start).count();
		double ontime_ms = std::chrono::duration_cast<std::chrono::milliseconds>(endOprf - startInputOn).count();

		std::cout << "\noffline= " << time1 << "ms\n";
		std::cout << "onlineCorrection= " << time2 << "ms\n";
		std::cout << "onlineOprf= " << time3 << "ms\n";
		std::cout << "total(online)= " << time_second << "(" << ontime_second << ") s\n";
		std::cout << "total(online)= " << time_ms << "(" << ontime_ms << ") ms\n\n";

		avg_total_time_ms += time_ms;
		avg_ontime_ms += ontime_ms;
#endif

#ifdef TEST
		recvChannel.recv((u8*)trans_pmSendMsg, sizeof(block)*(numBlkOPPRFs*numOTs));

		cout << "-- " << "0" << " " << "0" << ": " << trans_pmSendMsg[0] << "  " << trans_pmRecvMsg[0] << endl;


		for (u64 i = 0; i < numOPPRFs; ++i)
			for (u64 j = 0; j < numBlkPattern; ++j)
			{
				block revcBlock = trans_pmRecvMsg[i*numBlkPattern + j];
				block senderBlock = trans_pmSendMsg[i*numBlkPattern + j];

				if (eq(revcBlock, senderBlock))
					cout << "wildOn_trans-- " << i << " " << j << ": " << senderBlock << "  " << revcBlock << endl;

			}

		u64 i = 1, j = 2;
		block revcBlock = trans_pmRecvMsg[i*numBlkPattern + j];
		block senderBlock = trans_pmSendMsg[i*numBlkPattern + j];

		cout << "wildOn_trans-- " << i << " " << j << ": " << senderBlock << "  " << revcBlock << endl;

#endif // PRINT

		delete[] trans_pmSendMsg;
		delete[] trans_pmRecvMsg;
		delete[] blkPatternOff;
		delete[] blkPatternOn;
		delete[] blkWildOff;
		delete[] blkWildOn;
		for (u64 j = 0; j < 128; j += 1)
		{
			delete[] textShift[j];

		}


		recvChannel.close();
		ep0.stop(); ios.stop();
		}

	std::cout.imbue(std::locale(""));
	std::cout << "==========================\n ";
	std::cout << "|t|= " << numBlkText * 128 << "\t |p|=" << numBlkPattern * 128 << "\n";
	std::cout << "total/trial: total(online)= " << setprecision(2) << fixed << avg_total_time_ms / numTrials / 1000
		<< " (" << setprecision(2) << fixed << avg_ontime_ms / numTrials / 1000 << ") s\n";
	std::cout << "total/trial: total(online)= " << avg_total_time_ms / numTrials
		<< " (" << avg_ontime_ms / numTrials << ") ms\n";

	fOutputs << "==========================\n ";
	fOutputs << "|t|= " << numBlkText * 128 << "\t |p|=" << numBlkPattern * 128 << "\n";
	fOutputs << "total/trial: total(online)= " << setprecision(2) << fixed << avg_total_time_ms / numTrials / 1000
		<< " (" << setprecision(2) << fixed << avg_ontime_ms / numTrials / 1000 << ") s\n";
	fOutputs << "total/trial: total(online)= " << avg_total_time_ms / numTrials
		<< " (" << avg_ontime_ms / numTrials << ") ms\n";


	fLatex  << setprecision(2) << fixed << avg_total_time_ms / numTrials / 1000
		<< " (" << setprecision(2) << fixed << avg_ontime_ms / numTrials / 1000 << ")  &";

	fOutputs.close();
	fLatex.close();
	}

void usage(const char* argv0)
{
		std::cout << "Error! Please use:" << std::endl;
		std::cout << "\t 1. For unit test: " << argv0 << " -t" << std::endl;
		std::cout << "\t\t (default: nn=2^14, mm=2^9)" << std::endl;
		std::cout << "\t 2. For simulation (2 terminal): " << std::endl;;
		std::cout << "\t\t Sender terminal: " << argv0 << " -r 0 -nn 14 -mm 9" << std::endl;
		std::cout << "\t\t Receiver terminal: " << argv0 << "-r 1 -nn 14 -mm 9" << std::endl;
		std::cout << "\t\t (nn=2^16, mm=2^8)" << std::endl;
}


int main(int argc, char** argv)
{

	PRNG prng0(_mm_set_epi32(4253465, 3434565, 234435, 23987045));
	PRNG prng1(_mm_set_epi32(4253233465, 334565, 0, 235));
	
#if 0 //for debug
	{
#pragma region SetupInput
		u64 numBlkText = mNumBlkText; //
		u64 numBlkPattern = mNumBlkPattern; //32
		u64 numOTs = numBlkPattern * 128;
		u64 numBitText = numBlkText * 128;
		int numOPPRFs = numBitText - numOTs + 1;


		if (blk_step > numBlkPattern)
			blk_step = numBlkPattern;

		if (stepSize > numOPPRFs)
			stepSize = numOPPRFs;


		if (stepDeltaSize >numOTs)
			stepDeltaSize = numBlkPattern * 128;


		//Input
		BitVector  choicesWildOn(numOTs), choicesPatternOn(numOTs), choicesTextOn(numBitText);

		choicesWildOn.randomize(prng0);
		choicesPatternOn.randomize(prng0);
		choicesTextOn.randomize(prng0);

		std::array<int, 2> posMatches = { 0,numBitText - numOTs };


		for (u64 i = 0; i < numOTs; i++)
		{
			choicesTextOn[posMatches[0] + i] = choicesPatternOn[i];
			choicesTextOn[posMatches[1] + i] = choicesPatternOn[i];
		}

		for (u64 i = 0; i < 20; i++)
		{
			choicesTextOn[posMatches[1] + i] = choicesPatternOn[i] ^ 1;
			choicesWildOn[i] = 1;
		}

		cout << "|t|= " << numBitText << "\t |p|=" << numOTs << "\n";
		cout << "expected matching pos: " << posMatches[0] << ",\t" << posMatches[1] << "\n";

#pragma endregion

		std::thread thrd = std::thread([&]() {
			pmSender(choicesTextOn, numBlkText, numBlkPattern);
		});

		pmRecv(choicesPatternOn, choicesWildOn, numBlkText, numBlkPattern);

		thrd.join();

	}
	cin.ignore();
	return 0;
#endif

	if (argc == 2 && argv[1][0] == '-' && argv[1][1] == 't') {
			{
#pragma region SetupInput
			u64 numBlkText = mNumBlkText; //
			u64 numBlkPattern = mNumBlkPattern; //32
			u64 numOTs = numBlkPattern * 128;
			u64 numBitText = numBlkText * 128;
			int numOPPRFs = numBitText - numOTs+1;


			if (blk_step > numBlkPattern)
				blk_step = numBlkPattern;

			if (stepSize > numOPPRFs)
				stepSize = numOPPRFs;


			if (stepDeltaSize >numOTs)
				stepDeltaSize = numBlkPattern * 128;


			//Input
			BitVector  choicesWildOn(numOTs), choicesPatternOn(numOTs), choicesTextOn(numBitText);

			choicesWildOn.randomize(prng0);
			choicesPatternOn.randomize(prng0);
			choicesTextOn.randomize(prng0);

			std::array<int, 2> posMatches = { 0,numBitText - numOTs };


			for (u64 i = 0; i < numOTs; i++)
			{
				choicesTextOn[posMatches[0] + i] = choicesPatternOn[i];
				choicesTextOn[posMatches[1] + i] = choicesPatternOn[i];
			}

			for (u64 i = 0; i < 20; i++)
			{
				choicesTextOn[posMatches[1] + i] = choicesPatternOn[i] ^ 1;
				choicesWildOn[i] = 1;
			}

			cout << "|t|= " << numBitText << "\t |p|=" << numOTs << "\n";
			cout << "expected matching pos: " << posMatches[0] << ",\t" << posMatches[1] << "\n";

#pragma endregion

				std::thread thrd = std::thread([&]() {
					pmSender(choicesTextOn, numBlkText, numBlkPattern);
				});

				pmRecv(choicesPatternOn, choicesWildOn, numBlkText, numBlkPattern);

				thrd.join();

			}
	}
	else if (argc == 7 && argv[1][0] == '-' && argv[1][1] == 'r' && atoi(argv[2]) == 0) {
		cout << "sender\n";

#pragma region SetupInput
		u64 numBlkText,numBlkPattern;

		if(argv[3][0] == '-' && argv[3][1] == 'n' && argv[3][2] == 'n')
			numBlkText = 1 << (atoi(argv[4]) - 7);
		else if (argv[3][0] == '-' && argv[3][1] == 'n')
			numBlkText = ((atoi(argv[4]) + 127) / 128);
		else
		{
			usage(argv[0]);
			return 0;
		}

		if (argv[5][0] == '-' && argv[5][1] == 'm'&& argv[5][2] == 'm')
			numBlkPattern = 1 << (atoi(argv[6]) - 7);			
		else if (argv[5][0] == '-' && argv[5][1] == 'm')
			numBlkPattern = ((atoi(argv[6]) + 127) / 128);
		else
		{
			usage(argv[0]);
			return 0;
		}


		u64 numOTs = numBlkPattern * 128;
		u64 numBitText = numBlkText * 128;
		int numOPPRFs = numBitText - numOTs+1;


		if (blk_step > numBlkPattern)
			blk_step = numBlkPattern;

		if (stepSize > numOPPRFs)
			stepSize = numOPPRFs;


		if (stepDeltaSize >numOTs)
			stepDeltaSize = numBlkPattern * 128;

		//Input
		BitVector  choicesWildOn(numOTs), choicesPatternOn(numOTs), choicesTextOn(numBitText);

		choicesWildOn.randomize(prng0);
		choicesPatternOn.randomize(prng0);
		choicesTextOn.randomize(prng0);

		std::array<int, 2> posMatches = { 0,numBitText - numOTs };


		for (u64 i = 0; i < numOTs; i++)
		{
			choicesTextOn[posMatches[0] + i] = choicesPatternOn[i];
			choicesTextOn[posMatches[1] + i] = choicesPatternOn[i];
		}

		for (u64 i = 0; i < 20; i++)
		{
			choicesTextOn[posMatches[1] + i] = choicesPatternOn[i] ^ 1;
			choicesWildOn[i] = 1;
		}

		cout << "|t|= " << numBitText << "\t |p|=" << numOTs << "\n";
		cout << "expected matching pos: " << posMatches[0] << ",\t" << posMatches[1] << "\n";
#pragma endregion

		pmSender(choicesTextOn, numBlkText, numBlkPattern);
	}
	else if (argc == 7 && argv[1][0] == '-' && argv[1][1] == 'r' && atoi(argv[2]) == 1) {
		cout << "receiver\n";

#pragma region SetupInput
		u64 numBlkText, numBlkPattern;

		if (argv[3][0] == '-' && argv[3][1] == 'n' && argv[3][2] == 'n')
			numBlkText = 1 << (atoi(argv[4]) - 7);
		else if (argv[3][0] == '-' && argv[3][1] == 'n')
			numBlkText = ((atoi(argv[4]) + 127) / 128);
		else
		{
			usage(argv[0]);
			return 0;
		}

		if (argv[5][0] == '-' && argv[5][1] == 'm'&& argv[5][2] == 'm')
			numBlkPattern = 1 << (atoi(argv[6]) - 7);
		else if (argv[5][0] == '-' && argv[5][1] == 'm')
			numBlkPattern = ((atoi(argv[6]) + 127) / 128);
		else
		{
			usage(argv[0]);
			return 0;
		}


		u64 numOTs = numBlkPattern * 128;
		u64 numBitText = numBlkText * 128;
		int numOPPRFs = numBitText - numOTs + 1;


		if (blk_step > numBlkPattern)
			blk_step = numBlkPattern;

		if (stepSize > numOPPRFs)
			stepSize = numOPPRFs;


		if (stepDeltaSize >numOTs)
			stepDeltaSize = numBlkPattern * 128;


		//Input
		BitVector  choicesWildOn(numOTs), choicesPatternOn(numOTs), choicesTextOn(numBitText);

		choicesWildOn.randomize(prng0);
		choicesPatternOn.randomize(prng0);
		choicesTextOn.randomize(prng0);

		std::array<int, 2> posMatches = { 0,numBitText - numOTs };


		for (u64 i = 0; i < numOTs; i++)
		{
			choicesTextOn[posMatches[0] + i] = choicesPatternOn[i];
			choicesTextOn[posMatches[1] + i] = choicesPatternOn[i];
		}

		for (u64 i = 0; i < 20; i++)
		{
			choicesTextOn[posMatches[1] + i] = choicesPatternOn[i] ^ 1;
			choicesWildOn[i] = 1;
		}

		cout << "|t|= " << numBitText << "\t |p|=" << numOTs << "\n";
		cout << "expected matching pos: " << posMatches[0] << ",\t" << posMatches[1] << "\n";
#pragma endregion

		pmRecv(choicesPatternOn, choicesWildOn, numBlkText, numBlkPattern);
	}
	else {
		usage(argv[0]);
	}
    return 0;
}
