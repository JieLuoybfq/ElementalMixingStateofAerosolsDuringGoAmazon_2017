function sillystringout = sillystring

%LOTS OF VIDEO GAME SPOILERS AHEAD =P

sillystring = {...
    'Counting licks to center of tootsie pop.',... %tootsie pop commercial tag-line
    'Calculating exact value of pi.',... %made up
    'Determining quantum state of cat.',... %made up
    'Waking up the hamster.',... %made up
    'Testing your patience.',... %made up
    'This is still faster than you could do it.',... %made up
    'Loading humorous message.',... %made up
    'Insert quater to continue.',... %made up
    'Are we there yet?',... %made up
    'Help, im being forced to write these stupid lines',... %made up
    'Brewing coffee',... %made up
    'Dividing by 0',... %made up
    'If you had a better computer you wouldnt see this',... %made up
    'Loading virus and deleting important files',... %made up
    'Give me a second, I was sleeping',... %made up
    'The cake is a lie!',... %Rat Man scribbling on hidden wall in Portal 1
    'Maybe Black Mesa will help you. That was a joke, Ha Ha, fat chance',... %GlaDoS from Portal end credits song "Still Alive"
    'When life gives you lemons, don''t make lemonade. Make life take the lemons back!',... %Cave Johnson from Portal 2
    'Does this unit have a soul?',... %Legion from Mass Effect 3
    'Rise and shine, Mr. Freeman. Rise and shine.',... %The G-Man from Half-Life
    'Prepare for unforseen consequences',... %the G-Man/Alyx Vance from Half-Life 2
	'The right man in the wrong place can make all the difference in the world',... %The G-Man from Half-Life
    'Making a Jill sandwich...',... %Barry from Resident Evil 1
    'War... war never changes.',... %Intro to Fallout 1/2/3/NV/4 voiced by Ron Perlman
    'Your princess is in another castle',... %Toad from Super Mario Bros
    'Hey, listen!',... %Na'vi from Legend of Zelda Ocarina of Time
    'Stay a while, and listen.',... %Cain from Diablo 1/2/3
    'You have died of dysentery',... %Death message from The Oregon Trail
    'Would you kindly...',... %Atlas/Fontaine from Bioshock 1
    'A man chooses; a slave obeys.',... %Andrew Ryan from Bioshock 1
    'Is a man not entitled to the sweat of his brow?',... %Andrew Ryan from Bioshock 1
    'Gentlemen, welcome to Dubai',... %From Spec Ops: The Line
    'What is a man but a miserable little pile of secrets.',... %Dracula from Castlevania Symphony of the Night
    'What is a drop of rain compared to the storm?',... %The Many from System Shock 2
    'What is a thought compared to a mind?',... %The Many from System Shock 2
    'You must construct additional pylons',... %Protoss voice from Starcraft 1/2
    'Its dangerous to go alone! Take this.',... %Old Man from Legend of Zelda
    'Seeing this code run fills you with determination',... %Reference to save messages in Undertale
    'Fear the old blood',... %Warning given in Bloodborne by Master Willem to Laurence (Byrgenwerth's adage)
    'Seek Paleblood to transcend the hunt',... %Bloodborne random quote given near the beginning of the game
    'We are born by the blood, made men by the blood, undone by the blood',... % Master Willem from Bloodbourne
    'Soul of the mind, key to life''s ether. Soul of the lost, withdrawn from its vessel. Let strength be granted, so the world might be mended',... %Maiden in black from Demon's Souls
    'Lok''tar ogar',... %Orc saying from Warcraft series
    'The last metroid was in captivity, the galaxy is at peace',... %intro paragraph to Super Metroid
    'You cant ride your bike here',... %Prof. Oak in Pokemon games when trying to ride a bike indoors
    'In my restless dreams, I see that town... Silent Hill.  You promised you''d take me there... but you never did',... %Letter from Mary to James in Silent Hill 2
    'I look like Mary, don''t I? You loved her, right? Or maybe you hated her!',... %Maria talking to James about his dead wife from Silent Hill 2
    'You see it too? For me, it''s always like this...',... %Angela talking to James about the staircase burning around her from Silent Hill 2
    'This cannot continue. This cannot continue.',... %Machines just before combining to make/birth Adam from Nier Automata
    'Everything that lives is designed to end.  We are trapped in a never-ending spiral of life and death.  However, life is all about the struggle within this cycle.',... %Pod 153 ending E from Nier Automata
    'Bring us the girl and wipe away the debt',...}; %Bioshock Infinite

sillyidx = randi(length(sillystring));

sillystringout = sillystring{sillyidx};