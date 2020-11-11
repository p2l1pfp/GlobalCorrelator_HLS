library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use work.regionizer_data.all;

entity full_regionizer_mux_stream is
    generic(
            MU_ETA_CENTER : integer 
    );
    port(
            ap_clk : IN STD_LOGIC;
            ap_rst : IN STD_LOGIC;
            ap_start : IN STD_LOGIC;
            ap_done : OUT STD_LOGIC;
            ap_idle : OUT STD_LOGIC;
            ap_ready : OUT STD_LOGIC;
            tracks_start    : IN STD_LOGIC;
            tracks_newevent : IN STD_LOGIC;
            tracks_in_0_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_0_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_1_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_1_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_2_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_2_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_3_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_3_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_4_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_4_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_5_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_5_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_6_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_6_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_7_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_7_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_8_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            tracks_in_8_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_start    : IN STD_LOGIC;
            calo_newevent : IN STD_LOGIC;
            calo_in_0_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_0_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_0_2_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_0_3_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_1_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_1_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_1_2_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_1_3_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_2_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_2_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_2_2_V : IN STD_LOGIC_VECTOR (63 downto 0);
            calo_in_2_3_V : IN STD_LOGIC_VECTOR (63 downto 0);
            mu_start    : IN STD_LOGIC;
            mu_newevent : IN STD_LOGIC;
            mu_in_0_V : IN STD_LOGIC_VECTOR (63 downto 0);
            mu_in_1_V : IN STD_LOGIC_VECTOR (63 downto 0);
            --
            tracks_out : OUT w64s(NTKSTREAM-1 downto 0);
            calo_out   : OUT w64s(NCALOSTREAM-1 downto 0);
            mu_out     : OUT w64s(NMUSTREAM-1   downto 0);
            newevent_out : OUT STD_LOGIC
    );
end full_regionizer_mux_stream;

architecture Behavioral of full_regionizer_mux_stream is

    signal tracks_regionized:        w64s(NPFREGIONS-1 downto 0);
    signal tracks_regionized_valid:  std_logic_vector(NPFREGIONS-1 downto 0) := (others => '0');
    signal tracks_regionized_roll:   std_logic := '0';

    signal tracks_delayed:        w64s(TKDELAY*NPFREGIONS-1 downto 0);
    signal tracks_delayed_valid:  std_logic_vector(TKDELAY*NPFREGIONS-1 downto 0) := (others => '0');
    signal tracks_delayed_roll:   std_logic_vector(TKDELAY-1 downto 0);

    signal calo_regionized:        w64s(NPFREGIONS-1 downto 0);
    signal calo_regionized_valid:  std_logic_vector(NPFREGIONS-1 downto 0) := (others => '0');
    signal calo_regionized_roll:   std_logic := '0';

    signal mu_regionized:        w64s(NPFREGIONS-1 downto 0);
    signal mu_regionized_valid:  std_logic_vector(NPFREGIONS-1 downto 0) := (others => '0');
    signal mu_regionized_roll:   std_logic := '0';

    signal mu_delayed:        w64s(MUDELAY*NPFREGIONS-1 downto 0);
    signal mu_delayed_valid:  std_logic_vector(MUDELAY*NPFREGIONS-1 downto 0) := (others => '0');
    signal mu_delayed_roll:   std_logic_vector(MUDELAY-1 downto 0);

    signal tracks_sorted:        particles(NTKSORTED*NPFREGIONS-1 downto 0);
    signal tracks_sorted_valid:  std_logic_vector(NTKSORTED*NPFREGIONS-1 downto 0) := (others => '0');
    signal tracks_sorted_roll:   std_logic_vector(NPFREGIONS-1 downto 0) := (others => '0');

    signal calo_sorted:        particles(NCALOSORTED*NPFREGIONS-1 downto 0);
    signal calo_sorted_valid:  std_logic_vector(NCALOSORTED*NPFREGIONS-1 downto 0) := (others => '0');
    signal calo_sorted_roll:   std_logic_vector(NPFREGIONS-1 downto 0) := (others => '0');

    signal mu_sorted:        particles(NMUSORTED*NPFREGIONS-1 downto 0);
    signal mu_sorted_valid:  std_logic_vector(NMUSORTED*NPFREGIONS-1 downto 0) := (others => '0');
    signal mu_sorted_roll:   std_logic_vector(NPFREGIONS-1 downto 0) := (others => '0');

    signal tracks_mux :        particles(NTKSTREAM-1 downto 0);
    signal tracks_mux_valid :  std_logic_vector(NTKSTREAM-1 downto 0) := (others => '0');
    signal tracks_mux_roll :   std_logic := '0';

    signal calo_mux :        particles(NCALOSTREAM-1 downto 0);
    signal calo_mux_valid :  std_logic_vector(NCALOSTREAM-1 downto 0) := (others => '0');
    signal calo_mux_roll :   std_logic := '0';

    signal mu_mux :        particles(NMUSTREAM-1 downto 0);
    signal mu_mux_valid :  std_logic_vector(NMUSTREAM-1 downto 0) := (others => '0');
    signal mu_mux_roll :   std_logic := '0';

begin

    tk_regionizer : entity work.tk_regionizer 
                port map(ap_clk => ap_clk, ap_rst => ap_rst,
                             ap_start => tracks_start,
                             newevent => tracks_newevent,
                             tracks_in_0_0_V => tracks_in_0_0_V,
                             tracks_in_0_1_V => tracks_in_0_1_V,
                             tracks_in_1_0_V => tracks_in_1_0_V,
                             tracks_in_1_1_V => tracks_in_1_1_V,
                             tracks_in_2_0_V => tracks_in_2_0_V,
                             tracks_in_2_1_V => tracks_in_2_1_V,
                             tracks_in_3_0_V => tracks_in_3_0_V,
                             tracks_in_3_1_V => tracks_in_3_1_V,
                             tracks_in_4_0_V => tracks_in_4_0_V,
                             tracks_in_4_1_V => tracks_in_4_1_V,
                             tracks_in_5_0_V => tracks_in_5_0_V,
                             tracks_in_5_1_V => tracks_in_5_1_V,
                             tracks_in_6_0_V => tracks_in_6_0_V,
                             tracks_in_6_1_V => tracks_in_6_1_V,
                             tracks_in_7_0_V => tracks_in_7_0_V,
                             tracks_in_7_1_V => tracks_in_7_1_V,
                             tracks_in_8_0_V => tracks_in_8_0_V,
                             tracks_in_8_1_V => tracks_in_8_1_V,
                             tracks_out_0_V => tracks_regionized(0),
                             tracks_out_1_V => tracks_regionized(1),
                             tracks_out_2_V => tracks_regionized(2),
                             tracks_out_3_V => tracks_regionized(3),
                             tracks_out_4_V => tracks_regionized(4),
                             tracks_out_5_V => tracks_regionized(5),
                             tracks_out_6_V => tracks_regionized(6),
                             tracks_out_7_V => tracks_regionized(7),
                             tracks_out_8_V => tracks_regionized(8),
                             tracks_out_valid_0 => tracks_regionized_valid(0),
                             tracks_out_valid_1 => tracks_regionized_valid(1),
                             tracks_out_valid_2 => tracks_regionized_valid(2),
                             tracks_out_valid_3 => tracks_regionized_valid(3),
                             tracks_out_valid_4 => tracks_regionized_valid(4),
                             tracks_out_valid_5 => tracks_regionized_valid(5),
                             tracks_out_valid_6 => tracks_regionized_valid(6),
                             tracks_out_valid_7 => tracks_regionized_valid(7),
                             tracks_out_valid_8 => tracks_regionized_valid(8),
                             newevent_out => tracks_regionized_roll);


    calo_regionizer : entity work.calo_regionizer 
                port map(ap_clk => ap_clk, ap_rst => ap_rst,
                             ap_start => calo_start,
                             newevent => calo_newevent,
                             calo_in_0_0_V => calo_in_0_0_V,
                             calo_in_0_1_V => calo_in_0_1_V,
                             calo_in_0_2_V => calo_in_0_2_V,
                             calo_in_0_3_V => calo_in_0_3_V,
                             calo_in_1_0_V => calo_in_1_0_V,
                             calo_in_1_1_V => calo_in_1_1_V,
                             calo_in_1_2_V => calo_in_1_2_V,
                             calo_in_1_3_V => calo_in_1_3_V,
                             calo_in_2_0_V => calo_in_2_0_V,
                             calo_in_2_1_V => calo_in_2_1_V,
                             calo_in_2_2_V => calo_in_2_2_V,
                             calo_in_2_3_V => calo_in_2_3_V,
                             calo_out_0_V => calo_regionized(0),
                             calo_out_1_V => calo_regionized(1),
                             calo_out_2_V => calo_regionized(2),
                             calo_out_3_V => calo_regionized(3),
                             calo_out_4_V => calo_regionized(4),
                             calo_out_5_V => calo_regionized(5),
                             calo_out_6_V => calo_regionized(6),
                             calo_out_7_V => calo_regionized(7),
                             calo_out_8_V => calo_regionized(8),
                             calo_out_valid_0 => calo_regionized_valid(0),
                             calo_out_valid_1 => calo_regionized_valid(1),
                             calo_out_valid_2 => calo_regionized_valid(2),
                             calo_out_valid_3 => calo_regionized_valid(3),
                             calo_out_valid_4 => calo_regionized_valid(4),
                             calo_out_valid_5 => calo_regionized_valid(5),
                             calo_out_valid_6 => calo_regionized_valid(6),
                             calo_out_valid_7 => calo_regionized_valid(7),
                             calo_out_valid_8 => calo_regionized_valid(8),
                             newevent_out => calo_regionized_roll);

    mu_regionizer : entity work.mu_regionizer 
                generic map(ETA_CENTER => MU_ETA_CENTER)
                port map(ap_clk => ap_clk, ap_rst => ap_rst,
                             ap_start => mu_start,
                             newevent => mu_newevent,
                             mu_in_0_V => mu_in_0_V,
                             mu_in_1_V => mu_in_1_V,
                             mu_out_0_V => mu_regionized(0),
                             mu_out_1_V => mu_regionized(1),
                             mu_out_2_V => mu_regionized(2),
                             mu_out_3_V => mu_regionized(3),
                             mu_out_4_V => mu_regionized(4),
                             mu_out_5_V => mu_regionized(5),
                             mu_out_6_V => mu_regionized(6),
                             mu_out_7_V => mu_regionized(7),
                             mu_out_8_V => mu_regionized(8),
                             mu_out_valid_0 => mu_regionized_valid(0),
                             mu_out_valid_1 => mu_regionized_valid(1),
                             mu_out_valid_2 => mu_regionized_valid(2),
                             mu_out_valid_3 => mu_regionized_valid(3),
                             mu_out_valid_4 => mu_regionized_valid(4),
                             mu_out_valid_5 => mu_regionized_valid(5),
                             mu_out_valid_6 => mu_regionized_valid(6),
                             mu_out_valid_7 => mu_regionized_valid(7),
                             mu_out_valid_8 => mu_regionized_valid(8),
                             newevent_out => mu_regionized_roll);

    tk_delay : process(ap_clk)
            variable istart, iend: natural;
        begin
            if rising_edge(ap_clk) then
                istart := (TKDELAY-1)*NPFREGIONS; iend := TKDELAY*NPFREGIONS;
                tracks_delayed(      iend-1 downto istart) <= tracks_regionized(      NPFREGIONS-1 downto 0);
                tracks_delayed_valid(iend-1 downto istart) <= tracks_regionized_valid(NPFREGIONS-1 downto 0);
                tracks_delayed_roll(TKDELAY-1)             <= tracks_regionized_roll;
                if TKDELAY > 1 then
                    tracks_delayed(      istart-1 downto 0) <= tracks_delayed(      iend-1 downto NPFREGIONS);
                    tracks_delayed_valid(istart-1 downto 0) <= tracks_delayed_valid(iend-1 downto NPFREGIONS);
                    tracks_delayed_roll(TKDELAY-2 downto 0) <= tracks_delayed_roll(TKDELAY-1 downto 1);
                end if;
            end if;
       end process tk_delay;

   mu_delay : process(ap_clk)
            variable istart, iend: natural;
        begin
            if rising_edge(ap_clk) then
                istart := (MUDELAY-1)*NPFREGIONS; iend := MUDELAY*NPFREGIONS;
                mu_delayed(      iend-1 downto istart) <= mu_regionized(      NPFREGIONS-1 downto 0);
                mu_delayed_valid(iend-1 downto istart) <= mu_regionized_valid(NPFREGIONS-1 downto 0);
                mu_delayed_roll(MUDELAY-1)             <= mu_regionized_roll;
                if MUDELAY > 1 then
                    mu_delayed(      istart-1 downto 0) <= mu_delayed(      iend-1 downto NPFREGIONS);
                    mu_delayed_valid(istart-1 downto 0) <= mu_delayed_valid(iend-1 downto NPFREGIONS);
                    mu_delayed_roll(MUDELAY-2 downto 0) <= mu_delayed_roll(MUDELAY-1 downto 1);
                end if;
            end if;
       end process mu_delay;

    gen_sorters: for isort in NPFREGIONS-1 downto 0 generate
        tk_sorter : entity work.stream_sort
                            generic map(NITEMS => NTKSORTED)
                            port map(ap_clk => ap_clk,
                                d_in => w64_to_particle(tracks_delayed(isort)),
                                valid_in => tracks_delayed_valid(isort),
                                roll => tracks_delayed_roll(0),
                                d_out => tracks_sorted((isort+1)*NTKSORTED-1 downto isort*NTKSORTED),
                                valid_out => tracks_sorted_valid((isort+1)*NTKSORTED-1 downto isort*NTKSORTED),
                                roll_out => tracks_sorted_roll(isort)
                            );
        calo_sorter : entity work.stream_sort
                            generic map(NITEMS => NCALOSORTED)
                            port map(ap_clk => ap_clk,
                                d_in => w64_to_particle(calo_regionized(isort)),
                                valid_in => calo_regionized_valid(isort),
                                roll => calo_regionized_roll,
                                d_out => calo_sorted((isort+1)*NCALOSORTED-1 downto isort*NCALOSORTED),
                                valid_out => calo_sorted_valid((isort+1)*NCALOSORTED-1 downto isort*NCALOSORTED),
                                roll_out => calo_sorted_roll(isort)
                            );
        mu_sorter : entity work.stream_sort
                            generic map(NITEMS => NMUSORTED)
                            port map(ap_clk => ap_clk,
                                d_in => w64_to_particle(mu_delayed(isort)),
                                valid_in => mu_delayed_valid(isort),
                                roll => mu_delayed_roll(0),
                                d_out => mu_sorted((isort+1)*NMUSORTED-1 downto isort*NMUSORTED),
                                valid_out => mu_sorted_valid((isort+1)*NMUSORTED-1 downto isort*NMUSORTED),
                                roll_out => mu_sorted_roll(isort)
                            );
        end generate gen_sorters;

    tk_muxer: entity work.region_mux_stream
                            generic map(NREGIONS => NPFREGIONS, 
                                        NITEMS   => NTKSORTED,
                                        NSTREAM  => NTKSTREAM,
                                        OUTII    => PFII240)
                            port map(ap_clk => ap_clk,
                                roll => tracks_sorted_roll(0),
                                d_in => tracks_sorted,
                                valid_in => tracks_sorted_valid,
                                d_out => tracks_mux,
                                valid_out => tracks_mux_valid,
                                roll_out => tracks_mux_roll);

    calo_muxer: entity work.region_mux_stream
                            generic map(NREGIONS => NPFREGIONS, 
                                        NITEMS   => NCALOSORTED,
                                        NSTREAM  => NCALOSTREAM,
                                        OUTII    => PFII240)
                            port map(ap_clk => ap_clk,
                                roll => calo_sorted_roll(0),
                                d_in => calo_sorted,
                                valid_in => calo_sorted_valid,
                                d_out => calo_mux,
                                valid_out => calo_mux_valid,
                                roll_out => calo_mux_roll);

    mu_muxer: entity work.region_mux_stream
                            generic map(NREGIONS => NPFREGIONS, 
                                        NITEMS   => NMUSORTED,
                                        NSTREAM  => NMUSTREAM,
                                        OUTII    => PFII240)
                            port map(ap_clk => ap_clk,
                                roll => mu_sorted_roll(0),
                                d_in => mu_sorted,
                                valid_in => mu_sorted_valid,
                                d_out => mu_mux,
                                valid_out => mu_mux_valid,
                                roll_out => mu_mux_roll);

    format: process(ap_clk)
        begin
            if rising_edge(ap_clk) then
                for i in 0 to NTKSTREAM-1 loop
                    if tracks_mux_valid(i) = '1' then
                        tracks_out(i) <= particle_to_w64(tracks_mux(i));
                    else
                        tracks_out(i) <= (others => '0');
                    end if;
                end loop;
                for i in 0 to NCALOSTREAM-1 loop
                    if calo_mux_valid(i) = '1' then
                        calo_out(i) <= particle_to_w64(calo_mux(i));
                    else
                        calo_out(i) <= (others => '0');
                    end if;
                end loop;
                for i in 0 to NMUSTREAM-1 loop
                    if mu_mux_valid(i) = '1' then
                        mu_out(i) <= particle_to_w64(mu_mux(i));
                    else
                        mu_out(i) <= (others => '0');
                    end if;
                end loop;
                newevent_out <= tracks_mux_roll;
            end if;
        end process format;


end Behavioral;
